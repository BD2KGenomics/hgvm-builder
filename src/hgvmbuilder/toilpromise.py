#hgvm-builder toilpromise.py: Promise-like programming for Toil

import os
import os.path
import logging
import dill
import itertools

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
        # handlers will run as follow-ons. It can't be filled in until all our
        # handlers are set, because it pickles the handlers. This job will
        # return the promise's success result and its error, as a pair.
        self.promise_job = None
        
        # We might get an extra return value unpacking job
        self.unpack_resolve_job = None
        
        # This is the list of Promises that need to know about our Toil job and
        # so need to be start()'ed after we start() and make our job.
        self.dependents = []
        
        # We hold on to a Toil job to run as a child of, so we can schedule
        # ourselves after our handlers are added.
        self.parent_job = None
        
        # We could instead run as a follow-on of a Toil job.
        self.prev_job = None
        
        # We could instead have a Toil job we schedule ourselves after to
        # resolve with the return value.
        self.wrapped_job = None
        
        # We could instead be waiting for a single promise, which needs to be
        # .start()ed before us.
        self.prev_promise = None
        
        # We could instead be waiting for a bunch of promises, all of which need
        # to be .start()ed before us.
        self.wait_for = None
        
        # These are promises that are waiting on us, and which need to be
        # start()'ed when we are start()'ed so they can see if all their
        # dependencies have started yet.
        self.waiting_on_us = []
        
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
        
        # When this promise's executor calls its resolve, then we can run the
        # child's then handler.
        
        # Save the handler function
        child.set_then_handler(handler)
        
        # The child needs to be a follow-on of our job, when we get a job
        child.set_prev_promise(self)
        
        # And it's a dependent of us, so when we start it will start
        self.add_dependent(child)
        
        # Assuming we were started, start the child.
        child.start()
        
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
        
    def set_wrapped(self, job):
        """
        Wrap the given Toil job. Resolve when it finishes with its return value.
        We should not have an executor set.
        """
        
        # Make sure we have no real executor
        assert(self.executor_dill is None)
        
        self.wrapped_job = job
        
        
    def set_parent_job(self, parent_job):
        """
        Make the Promise be a child of the given parent Toil job, where it will
        run its executor or resolve to its set value.
        
        """
        
        assert(self.executor_dill is not None or self.result is not None)
        
        self.parent_job = parent_job
        
    def set_prev_job(self, prev_job):
        """
        Make the Promise be a follow-on of the given predecessor Toil job, where
        it will run its executor or resolve to its set value.
        
        """
        
        assert(self.executor_dill is not None or self.result is not None)
        
        self.prev_job = prev_job
        
    def set_prev_promise(self, prev_promise):
        """
        Make the Promise be a follow-on of the given predecessor promise, where
        it will run its then handler.
        
        """
        
        # We must have a then handler
        assert(self.then_dill is not None)
        
        self.prev_promise = prev_promise
        
    def set_wait_for(self, promises):
        """
        Set this Promise to wait for all the given promises and then resolve
        with their results.
        """
        
        self.wait_for = promises
        
    def add_dependent(self, promise):
        """
        Remember to start() the given promise after we start() and make our Toil
        job to depend on.
        """
        
        self.dependents.append(promise)
        
    def start(self):
        """
        We're done adding handlers, so pickle the promise and schedule it in
        Toil. Also make sure to schedule anything that depends on it.
        
        All the promises need to be start()ed in the Toil job that defined the
        whole promise graph, because they need to create their jobs so they can
        be linked up in the dependency graph.
        
        """
        
        # TODO: This really ought to be polymorphism based on type instead of
        # dispatch based on filled fields.
        
        if self.promise_job is None:
            # Haven't made a Toil job yet, so try it.
        
            if self.parent_job is not None and self.result is not None:
                # We want a promise that always resolves, as a child job.
                self.promise_job = self.parent_job.addChildJobFn(
                    promise_resolve_job, self)
                    
            elif self.prev_job is not None and self.result is not None:
                # We want a promise that always resolves, as a followon job.
                self.promise_job = self.prev_job.addFollowOnJobFn(
                    promise_resolve_job, self)
                    
            elif self.parent_job is not None and self.executor_dill is not None:
                # We want to run this executor closure under this job as a child
                self.promise_job = self.parent_job.addChildJobFn(
                    promise_executor_job, self)
                
            elif self.prev_job is not None and self.executor_dill is not None:
                # We want to run this executor closure under this job as a
                # follow-on
                self.promise_job = self.prev_job.addFollowOnJobFn(
                    promise_executor_job, self)
                    
            elif self.prev_promise is not None and self.then_dill is not None:
                # We are handling a .then() promise that depends on one parent.
                # We assume the parent is started, and add as a follow-on.
                self.promise_job = \
                    self.prev_promise.promise_job.addFollowOnJobFn(
                    promise_then_job, self, self.prev_promise.promise_job.rv())
                    
            elif self.wrapped_job is not None and self.executor_dill is None:
                # We want to wrap this Toil job
                self.promise_job = self.wrapped_job.addFollowOnJobFn(
                    promise_wrapper_job, self, self.wrapped_job.rv())
                    
            elif self.wait_for is not None and self.executor_dill is None:
                # We want to run after all these started promises and resolve with
                # their results.
                
                
                # Have all the promises we need to wait for been start()'ed?
                all_started = True
                
                # Make a list of all the Toil RVs to wait for. Each is a pair of
                # result and error from the promises we are waiting on.
                rv_list = []
                
                for waited in self.wait_for:
                    if waited.promise_job is not None:
                        # This dependency has start()'ed
                        rv_list.append(waited.promise_job.rv())
                    else:
                        all_started = False
                        
                if all_started:
                    # We can actually schedule!
                    
                    # Make a job to resolve us with all the successes, or reject
                    # us with the failures
                    self.promise_job = Job.wrapJobFn(promise_all_job, self,
                        rv_list)
                    
                    # Schedule it after all the promises it is supposed to wait
                    # for.
                    for waited in self.wait_for:
                        waited.promise_job.addFollowOn(self.promise_job)
                
            else:
                raise RuntimeError("Invalid promise configuration!")
          
        if self.promise_job is not None:
            # Our dependencies are all ready, so we are too!
            for dependent in self.dependents:
                # Start the dependents, which may have jobs created for all their
                # dependsencies now.
                dependent.start()
                
        # TODO: everything should be started when it's made, since then handlers
        # now live in their own promises. We can basically cut this function!
        assert(self.promise_job is not None)
            
        
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
        After a promise has been .start()'ed, you can get the Toil return value
        .rv() for its resolve-value, reject-value pair.
        
        You can use this to convert from promises back to normal Toil style
        code.
        
        Should only be used in the parent job of the entire promise tree.
        
        """
        
        assert(self.promise_job is not None)
        
        return self.promise_job.rv()
        
    def unwrap_result(self):
        """
        After a promise has been .start()'ed, you can get the Toil return value
        .rv() for just its resolve value, or None if it rejects.
        
        You can use this to convert from promises back to normal Toil style
        code.
        
        Should only be used in the parent job of the entire promise tree.
        
        """
        
        assert(self.promise_job is not None)
        
        if self.unpack_resolve_job is None:
            # We need to add an extra follow-on to grab the resolve value (0th returned)
            self.unpack_resolve_job = self.promise_job.addFollowOnJobFn(
                index_job, self.promise_job.rv(), 0)
                
        # Once that job exists once, we can just grab its rv().
        return self.unpack_resolve_job.rv()
        
    @staticmethod
    def resolve(parent_job, result):
        """
        Produce a Promise that always resolves with the given result when run.
        
        The result cannot contain Toil .rv()s from children of the parent job,
        or you will get EOFErrors from Toil's unpickling code! Access .rv()s in
        Promises with wrap!
        
        Promise will need to be manually .start()'ed after handlers are added.
        
        """
        
        # Don't use a lambda to always resolve, because the thing we want to
        # resolve to may have .rv()s in it. Send it via normal Toil
        # serialization.
        
        promise = ToilPromise()
        promise.result = result
        promise.set_parent_job(parent_job)
        
        # Since the parent job exists, we can make our job.
        promise.start()
        
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
        promise.set_prev_job(prev_job)
        
        # Since the prev job exists, we can make our job.
        promise.start()
        
        return promise
            
    @staticmethod
    def all(promises):
        """
        Produce a Promise that resolves when all the given promises resolve. The
        promise resolves with all their results in a list.
        
        """
    
        # Make a Promise
        after_promise = ToilPromise()
    
        # We can't set it as a follow-on of anything until it has a job, and it
        # can't have a job until it gets started in start. So we just remember
        # the promises to wait for.
        after_promise.set_wait_for(promises)
        
        # But those promises won't have jobs to wait on until they start, so we
        # need to register to be started by them.
        for promise in promises:
            promise.add_dependent(after_promise)
        
        # Assuming all the promises we wanted to depend on have started, start
        # this one.
        after_promise.start()
        
        return after_promise
        
    @staticmethod
    def wrap(toil_job):
        """
        Produce a Promise that wraps the given normal Toil job, and resolves
        with its return value.
        
        Promise will need to be manually .start()'ed after handlers are added.
        
        """
        
        wrapper_promise = ToilPromise()
        wrapper_promise.set_wrapped(toil_job)
        
        # Since the wrapped job exists, we can make our job.
        wrapper_promise.start()
        
        return wrapper_promise
        
    @staticmethod
    def add_child(parent_job, executor):
        """
        Make a promise as a child of a Toil job, running the given executor.
        
        The executor is a function that calls its first argument with a result,
        or its second argument with an error.
        
        Promise will need to be manually .start()'ed after handlers are added.
        
        """
        
        child_promise = ToilPromise()
        child_promise.set_executor(executor)
        child_promise.set_parent_job(parent_job)
        
        # Since the parent job exists, we can make our job.
        child_promise.start()
        
        return child_promise
        
    @staticmethod
    def add_followon(prev_job, executor):
        """
        Make a promise as a follow-on of a Toil job, running the given executor.
        
        The executor is a function that calls its first argument with a result,
        or its second argument with an error.
        
        Promise will need to be manually .start()'ed after handlers are added.
        
        """
        
        child_promise = ToilPromise()
        child_promise.set_executor(executor)
        child_promise.set_prev_job(prev_job)
        
        # Since the prev job exists, we can make our job.
        child_promise.start()
        
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
    Make the given promise resolve with the list of the first elements in all
    the pairs, if the second elements are all none, or reject with a ist of the
    second elements, if any are not None.
    
    Makes ToilPromise.all() work.
    """
    
    # Looka t all the results
    results = []
    errs = []
    failed = False
    
    for result, err in result_pairs:
        # Put the results and errors in the right list
        results.append(result)
        errs.append(err)
        
        if err is not None:
            # If we have any errors we should fail
            failed = True

    if failed:
        # Fail the promise if anything in the list failed
        promise.handle_reject(job, errs)
    else:
        # Succeed otherwise
        promise.handle_resolve(job, results)
        
    # Return what the promise resolved/rejected with
    return (promise.result, promise.err)
    
def index_job(job, collection, item):
    """
    Grab the item at the given index and return it, as a Toil job.
    """
    
    return collection[item]




















