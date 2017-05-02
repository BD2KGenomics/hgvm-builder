#hgvm-builder toilpromise.py: Promise-like programming for Toil

import os
import os.path
import logging
import dill
import itertools
import sys
import traceback

from toil.job import Job
from toil.realtimeLogger import RealtimeLogger

Logger = logging.getLogger("toilpromise")

class ToilPromise(object):
    """
    Represents a Promise. A promise (theoretically) has an executor function
    that takes resolve and reject functions as arguments. In practice we have a
    bunch of fake executor functions.
    
    When the executor succeeds, it must call its resolve function with its
    result. When it fails, it must call its reject function with its error
    message.
    
    You can use .then(success_handler) to register a function to be called with
    the result when the promise succeeds.
    
    Handlers will be called in follow-on Toil jobs, in order to resolve their
    return values.
    
    TODO: We should have a base class and a bunch of polymorphic implementations
    that all implement something like .get_toil_job().
    
    """
    
    def __init__(self):
        """
        Make a new Promise that does nothing. Not to be used outside this class.
        
        """
        
        # These fields pickle
        
        # We have no executor yet
        self.executor_dill = None
        
        # We also have no then handler, which is run as a fake executor if we
        # have it
        self.then_dill = None
        
        # When the promise resolves locally, it fills in this result. Can be
        # filled in for a promise that always resolves.
        self.result = None
        
        # When it errors locally, it fills in this error.
        self.err = None
        
        # These fields don't pickle
        
        # This Toil job represents the actual execution of the promise. All the
        # handlers will run as follow-ons. This job will return the promise's
        # success result and its error, as a pair.
        self.promise_job = None
        
        # We might get an extra return value unpacking job
        self.unpack_resolve_job = None
        
    def __getstate__(self):
        """
        Return a tuple of state to be pickled.
        """
        
        return (self.then_dill, self.executor_dill, self.result, self.err)

    def __setstate__(self, state):
        """
        Restore state from a tuple.
        """
        
        (self.then_dill, self.executor_dill, self.result, self.err) = state
        
        
    def then(self, handler):
        """
        Schedule the given handler to be called with the resolve result.
        
        Returns a promise that resolves with the handler's return value.
        """
        
        # Make the new promise
        child = ToilPromise()
        # Save the handler function
        child.set_then_handler(handler)
        
        # Make the promise's job as a follow-on of ours, getting our return
        # value.
        child.promise_job = self.promise_job.addFollowOnJobFn(promise_then_job,
            child, self.promise_job.rv())
        
        return child
        
    def then_job_fn(self, job_fn, *args, **kwargs):
        """
        Pass the result of this ToilPromise to the given Toil job function,
        which takes a Toil job as a first argument.
        
        Additional arguments after the job function name will be passed before
        the promise's resolve value.
        
        If the promise rejects, the added Toil job may still be called, with a
        None argument.
        
        Makes sure to pickle the function, so it can be a lambda.
        
        """
        
        # Make the new promise
        child = ToilPromise()
        # Save the handler function
        child.set_then_handler(job_fn)
        
        # Stick the result, error pair from the last promise at the end of any
        # forwarded args, where the handler expects the result to end up.
        args = list(args) + [self.promise_job.rv()]
        
        
        # Make the promise's job as a follow-on of ours, getting our return
        # value.
        child.promise_job = self.promise_job.addFollowOnJobFn(
            promise_then_job_fn_job, child, *args, **kwargs)
        

        # Return the child, which looks like a then handler but is associated
        # with a different kind of Toil job.
        return child
        
    def set_executor(self, executor):
        """
        Set the given function to run in the promise. It will call its first
        argument with its result, or its second argument with an error.
        """
        
        # Pickle the function and save it
        self.executor_dill = dill.dumps(executor)
        
    def set_then_handler(self, then_handler):
        """
        Set the then handler for this promise. When the prev promise resolves,
        the then handler will be called with the result.
        
        """
        
        # Pickle the function and save it
        self.then_dill = dill.dumps(then_handler)
        
    def handle_resolve(self, job, result):
        """
        Handle promise resolution. Runs in the promise's assigned job.
        """
        
        Logger.info("Promise Resolved: {}".format(result))
        
        # Cache the result
        self.result = result
        
        # Toil already knows to run the Promises that depend on this result.
            
    def handle_reject(self, job, err):
        """
        Handle promise rejection.
        """
        
        Logger.error("Promise Rejected: {}".format(err))
        RealtimeLogger.error("Promise Rejected: {}".format(err))
        
        self.err = err
        
        # TODO: implement
        # Check if we have any reject handlers
        # If so call them
        # If not throw an error that stops the workflow
        raise err
        
    def unwrap(self):
        """
        You can get the Toil return value .rv() for its resolve-value, reject-
        value pair.
        
        You can use this to convert from promises back to normal Toil style
        code.
        
        Should only be used in the parent job of the entire promise tree.
        
        """
        
        assert(self.promise_job is not None)
        
        return self.promise_job.rv()

    def unwrap_result(self, *args):
        """
        You can get the Toil return value .rv() for just its resolve value, or
        None if it rejects.
        
        This is basically the same as .unwrap()(1), using Toil's built-in
        promise indexing, but with some extra checks.
        
        You can use this to convert from promises back to normal Toil style
        code, but if a converted promise rejects, it will kill the workflow.
        
        You can also pass a toil-style path (anything you would pass to .rv())
        and it will be used to pull stuff out of the result that got unwrapped.
        
        Should only be used in the parent job of the entire promise tree.
        
        """
        
        assert(self.promise_job is not None)
        
        if self.unpack_resolve_job is None:
            # We need to add an extra follow-on to grab the resolve value, or
            # fail if it rejected.
            self.unpack_resolve_job = self.promise_job.addFollowOnJobFn(
                get_resolve_value_job, self.promise_job.rv())
                
        # Once that job exists once, we can just grab its rv().
        return self.unpack_resolve_job.rv(*args)
        
    def get_last_job(self):
        """
        Grab the Toil job that you need to be a follow-on of in order to safely
        use the result of .unwrap() or .unwrap_result().
        
        You might want to use addFollowOnJobFn instead if you want Toil code to
        wait on a promise.
        
        """
        
        # Because unwrap_result uses another joib, we return that one.
        
        # Make sure it exists
        self.unwrap_result()
        # Return it
        return self.unpack_resolve_job
        
    def addFollowOnJobFn(self, *args, **kwargs):
        """
        Add a Toil job as a follow-on of the last job that can possibly belong
        to this promise. That Toil follow on will be able to read the result of
        .unwrap() or .unwrap_result().
        
        """
        
        return self.get_last_job().addFollowOnJobFn(*args, **kwargs)
        
    @staticmethod
    def resolve(parent_job, result):
        """
        Produce a Promise that always resolves with the given result when run.
        
        The result cannot contain Toil .rv()s from children of the parent job,
        or you will get EOFErrors from Toil's unpickling code! Access .rv()s in
        Promises with wrap!
        
        """
        
        # Don't use a lambda to always resolve, because the thing we want to
        # resolve to may have .rv()s in it. Send it via normal Toil
        # serialization.
        
        promise = ToilPromise()
        promise.result = result
        
        # Make the promise's job
        promise.promise_job = parent_job.addChildJobFn(promise_resolve_job,
            promise)
        
        return promise
        
    @staticmethod
    def resolve_followon(prev_job, result):
        """
        Produce a Promise that always resolves with the given result when run.
        The Promise will run as a follow-on of the given job.
        
        Using this on the current Toil job can make cycles in the Toil graph
        (somehow)!
        
        """
        
        promise = ToilPromise()
        promise.result = result
        
        # Make the promise's job
        promise.promise_job = prev_job.addFollowOnJobFn(promise_resolve_job,
            promise)
        
        return promise
            
    @staticmethod
    def all(promises):
        """
        Produce a Promise that resolves when all the given promises resolve.
        Works on either a list of promises, or a dict of promises by some kind
        of value key, and resolves to either a list of result values/reject
        errors, or a dict of result values/reject errors by key.
        
        """
    
        # Make a Promise
        after_promise = ToilPromise()
        
        if isinstance(promises, list):
            # Make a list of RVs
            rv_struct = [waited.promise_job.rv() for waited in promises]
            # And a set of unique jobs
            unique_jobs = {waited.promise_job for waited in promises}
        elif isinstance(promises, dict):
            # Make a dict of RVs
            rv_struct = {k: waited.promise_job.rv()
                for k, waited in promises.iteritems()}
            # And a set of unique jobs
            unique_jobs = {waited.promise_job for waited in promises.values()}
        else:
            raise RuntimeError("Unsupported data structure for waiting on")
            
        
         # Make a job to represent the promise, passing the data structure that
         # will get filled in with the (resolve value, reject error) pairs.
        after_promise.promise_job = Job.wrapJobFn(promise_all_job,
            after_promise, rv_struct)
            
        for waited in unique_jobs:
            # Run the promise's job after all the jobs it depends on
            waited.addFollowOn(after_promise.promise_job)
        
        return after_promise
        
    @staticmethod
    def wrap(toil_job):
        """
        Produce a Promise that wraps the given normal Toil job, and resolves
        with its return value.
        
        """
        
        wrapper_promise = ToilPromise()
        
        # Make the promise's job
        wrapper_promise.promise_job = toil_job.addFollowOnJobFn(
            promise_wrapper_job, wrapper_promise, toil_job.rv())
        
        return wrapper_promise
        
    @staticmethod
    def add_child(parent_job, executor):
        """
        Make a promise as a child of a Toil job, running the given executor.
        
        The executor is a function that calls its first argument with a result,
        or its second argument with an error.
        
        """
        
        child_promise = ToilPromise()
        child_promise.set_executor(executor)

        # Make the promise's job        
        child_promise.promise_job = parent_job.addChildJobFn(
            promise_executor_job, child_promise)
        
        return child_promise
        
    @staticmethod
    def add_followon(prev_job, executor):
        """
        Make a promise as a follow-on of a Toil job, running the given executor.
        
        The executor is a function that calls its first argument with a result,
        or its second argument with an error.
        
        """
        
        child_promise = ToilPromise()
        child_promise.set_executor(executor)

        # Make the promise's job        
        child_promise.promise_job = prev_job.addFollowOnJobFn(
            promise_executor_job, child_promise)
        
        return child_promise
        

def promise_resolve_job(job, promise):
    """
    Toil job that resolves a promise with its result and err already filled.
    """
    
    if promise.err is None:
        promise.handle_resolve(job, promise.result)
    else:
        promise.handle_reject(job, promise.err)
        
    return (promise.result, promise.err)

def promise_executor_job(job, promise):
    """
    Toil job that runs a promise with an executor. Executes the executor and
    rejects/resolves the promise.
    
    Returns the promise's success result and error, as a pair.
    
    """

    executor = dill.loads(promise.executor_dill)
    
    # Run the executor, and handle resolution/rejection, possibly scheduling
    # child jobs
    executor(lambda result: promise.handle_resolve(job, result),
        lambda err: promise.handle_reject(job, err))
        
    # Grab the cached result and return it
    return (promise.result, promise.err)
    
def promise_then_job(job, promise, prev_promise_returned):
    """
    Toil job that runs a promise created with a then handler instead of an
    executor.
    
    Takes the promise and the (resolve value, reject value) pair from the
    previous promise.
    
    Returns the promise's success result and error, as a pair.
    """
    
    then_handler = dill.loads(promise.then_dill)
    
    resolved, rejected = prev_promise_returned
    
    if rejected is None:
        # Actually run this child promise
        
        try:
            # Get the result from the then handler and resolve with it
            result = then_handler(resolved)
            promise.handle_resolve(job, result)
        except Exception as e:
            # Reject with an error if there is one
            Logger.error("".join(traceback.format_exception(*sys.exc_info())))
            promise.handle_reject(job, e)
            
    else:
        # Parent promise rejected so we should not run
        # Bubble up the error
        promise.handle_reject(job, rejected)
        
    return (promise.result, promise.err)
    
def promise_then_job_fn_job(job, promise, *args, **kwargs):
    """
    Toil job that runs a promise created with a then_job_fn handler.
    
    Takes the promise, and the arguments to forward along to the handler, the
    last of which is the (result, error) pair from the last promise which gets
    processed to just a result.
    
    Returns the promise's success result and error, as a pair.
    """
    
    
    # The pickled handler in this case takes a bunch of arguments: the Toil job,
    # and the success result from the last promise, and then any other arguments
    # or kwargs that the user wanted to pass along.
    then_handler = dill.loads(promise.then_dill)
    
    # Pull out the results from the last promise
    resolved, rejected = args[-1]
    args = args[:-1]
    
    if rejected is None:
        # Actually run this child promise
        
        # Stick the resolved value on
        args = list(args) + [resolved]
        
        try:
            # Get the result from the then handler and resolve with it
            result = then_handler(job, *args, **kwargs)
            promise.handle_resolve(job, result)
        except Exception as e:
            # Reject with an error if there is one
            Logger.error("".join(traceback.format_exception(*sys.exc_info())))
            promise.handle_reject(job, e)
            
    else:
        # Parent promise rejected so we should not run
        # Bubble up the error
        promise.handle_reject(job, rejected)
        
    return (promise.result, promise.err)
        
def promise_wrapper_job(job, promise, result):
    """
    Make the given promise resolve with the given result, which was a Toil-
    internal .rv() promise, but which is now real.
    
    Returns the result again, and None for the error, in a pair.
    """
    
    promise.handle_resolve(job, result)
    
    return (promise.result, promise.err)
    
def promise_all_job(job, promise, result_pairs):
    """
    Operates on a data structure (list or dict) with (resolve result, reject
    error) pairs as values.
    
    Make the given promise resolve with the struct of the first elements in all
    the pairs, if the second elements are all none, or reject with a struct of
    the second elements, if any are not None.
    
    Makes ToilPromise.all() work.
    """
    
    # Did anything we depend on fail?
    failed = False
    
    if isinstance(result_pairs, list):
        # Swizzle out into lists
        results = []
        errs = []
        
        for result, err in result_pairs:
            # Put the results and errors in the right list
            results.append(result)
            errs.append(err)
            
            if err is not None:
                # If we have any errors we should fail
                failed = True
                
    elif isinstance(result_pairs, dict):
        # Swizzle out into dicts
        results = {}
        errs = {}
        
        for key, (result, err) in result_pairs.iteritems():
            # Put the results and errors in the right list
            results[key] = result
            errs[key] = err
            
            if err is not None:
                # If we have any errors we should fail
                failed = True
    else:
        raise RuntimeError("Unsupported data structure!")

    if failed:
        # Fail the promise if anything in the list failed
        promise.handle_reject(job, errs)
    else:
        # Succeed otherwise
        promise.handle_resolve(job, results)
        
    # Return what the promise resolved/rejected with
    return (promise.result, promise.err)
    
def get_resolve_value_job(job, resolve_reject):
    """
    Grab the resolve value from a promise job's resolve, reject return pair.
    Fail if the promise rejected.
    """
    
    # Make sure the promise resolved
    assert(resolve_reject[1] is None)
    # Unpack the resolved value
    return resolve_reject[0]




















