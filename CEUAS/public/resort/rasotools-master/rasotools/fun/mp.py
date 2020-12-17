#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import subprocess
from time import time, sleep


class SubProcess:
    def __init__(self, icmd, shell=False, ref=None):
        self.proc = subprocess.Popen(icmd.split() if not shell else icmd, stdout=subprocess.PIPE, shell=shell)
        self.start = time()
        self.cmd = icmd
        self.ref = ref

    def check(self, timeout=None):
        if self.proc.poll() is None:
            # running
            if timeout is not None:
                # timeout
                if (time() - self.start) > timeout:
                    # over the limit
                    # print("Terminated", time() - self.start, self.cmd)
                    self.proc.terminate()
                    sleep(0.1)
                    return False

            return True
        else:
            return False

    def get(self):
        if self.proc.poll() is not None:
            if self.proc.poll() < 0:
                return "Error " + str(self.proc.poll()) + self.proc.stdout.read().decode() + '\n'  # add linefeed ?
            else:
                return self.proc.stdout.read().decode()
        return None


def batch_shell_command(command, files, max_processes=10, shell=False, interval=1, timeout=None, return_list=False):
    """ Run a shell command on multiple files an get results back

    Parameters
    ----------
        command : str
            Command string including FILE as replacent
        files : list
            list of files
        max_processes : int
            number of process to run
        shell : bool
            execute as a shell command (including pipes)
        interval : float
            sleep time
        timeout : float
            maximum allowed time for process to run (seconds)

    Returns
    -------
    list
        List of Results
    """
    processes = set()
    data = []
    indices = []
    n = len(files)
    i = 0
    m = 0
    debug = False

    while True:
        while len(processes) < max_processes and i < n:
            # starte neue
            name = files[i]
            icmd = command.replace('FILE', name)
            processes.add(SubProcess(icmd, shell=shell, ref=name))
            if debug: print("Launched", icmd)
            i += 1
            k = round((i / n) * 100, 1)
            if k != m:
                print("Status: {:04.1f} %".format(k))  # Progress report
                m = k

        sleep(interval)
        done = []
        for p in processes:
            if debug: print("Current", p.cmd, time() - p.start)

            if not p.check(timeout=timeout):
                data.append(p.get())
                done.append(p)
                indices.append(p.ref)

        processes.difference_update(done)
        if i >= n:
            break

    if debug: print("Collecting ...")
    while len(processes) > 0:
        done = []
        for p in processes:
            if not p.check(timeout=timeout):
                if debug: print("Done", p.cmd, time() - p.start)
                data.append(p.get())
                done.append(p)
                indices.append(p.ref)
        processes.difference_update(done)
        if len(processes) > 0:
            sleep(interval)

    if return_list:
        return data, indices
    return data


class MPProcess:
    def __init__(self, cmd, name, *args, **kwargs):
        self.cmd = cmd
        self.name = name
        self.arguments = args
        self.keywords = kwargs
        # Return stuff
        self.executed = False
        self.return_value = False
        self.status = False
        self.error = ''
        self.value = None

    def __repr__(self):
        out = "<Process %s >\n" % self.name
        out += "\n".join({"%s : %s " % (i, repr(j)) for i, j in self.__dict__.items() if i != 'name'})
        return out

    def execute(self):
        if not self.executed:
            self.value = None
            self.return_value = self.keywords.pop('return_value', False)  # removed

            try:
                self.value = self.cmd(*self.arguments, **self.keywords)
                self.status = True
                self.error = ""
            except Exception as e:
                self.status = False
                self.error = repr(e)

            self.executed = True
        return self


def worker(x):
    return x.execute()


def make_process_list(iterable, function, args=(), **kwargs):
    stuff = []
    for i in iterable:
        iargs = (i,) + args
        stuff.append(MPProcess(function, i, *iargs, **kwargs))
    return stuff


def execute_process(iterable, npp=10, return_value=False, loop=False, progress=True, summary=True, extend=False):
    try:
        import tqdm
    except ImportError:
        progress = False
        print("No tqdm module available, no progress update")
    try:
        from multiprocessing import Pool
    except ImportError:
        loop = True
        print("No multiprocessing available, using loop")

    if not isinstance(iterable[0], MPProcess):
        raise ValueError("Requires a MPProcess object")

    res = []
    if loop:
        if progress:
            for p in tqdm.tqdm(iterable, total=len(iterable)):
                res.append(p.execute())
        else:
            for p in iterable:
                res.append(p.execute())
    else:
        try:
            if progress:
                with Pool(npp) as p:
                    res = list(tqdm.tqdm(p.imap(worker, iterable), total=len(iterable)))
            else:
                with Pool(npp) as p:
                    res = p.imap(worker, iterable)

        except Exception as e:
            print("Error in MP", repr(e))

        else:
            p.terminate()

    # todo what if not all have been executed ? / Execute others?
    # p.execute() if not self.executed ?
    results = []
    ok = []
    errors = []
    for i in res:
        if i.status and i.executed:
            ok.append(i.name)
            if extend:
                results.extend(i.value)
            else:
                results.append(i.value)
        else:
            errors.append(i.name)
            print(i.error)
    if summary:
        print("MP [%d] %s Status: %4.0f %% (# Failure: %d)" % (
            npp, 'LP' if loop else 'MP', 100 * len(iterable) / len(ok), len(errors)))

    if return_value:
        return errors, list(results)


if __name__ == "__main__":
    import numpy as np

    #
    # Example
    #
    ivals = np.random.rand(10000)
    todo = make_process_list(ivals, np.cos, return_value=True)

    print("Executing Cos on 10000 random numbers")
    td1 = time()
    print(np.cos(ivals))
    print("Time (np): ", time() - td1)
    td1 = time()
    print(execute_process(todo, 10, return_value=True, loop=True))
    print("Time (loop): ", time() - td1)
    td1 = time()
    print(execute_process(todo, 10, return_value=True))
    print("Time (pool): ", time() - td1)


    def test_function(x):
        return np.sum(x)


    print("Executing sum of [0:n] n=1000000")
    todo = make_process_list([np.random.rand(n) for n in range(1, 1000000)], test_function, return_value=True)
    td1 = time()
    execute_process(todo, 10, return_value=False, loop=True)
    print("Time (loop): ", time() - td1)
    td1 = time()
    execute_process(todo, 10, return_value=False)
    print("Time (pool): ", time() - td1)


    #
    #
    #
    # Calculate Fibunacci numbers
    #
    #
    #
    #

    def fibonacci(n):
        if n < 0:
            print("Incorrect input")
            # First Fibonacci number is 0
        elif n == 1:
            return 0
        # Second Fibonacci number is 1
        elif n == 2:
            return 1
        else:
            return fibonacci(n - 1) + fibonacci(n - 2)


    print("Fibunacci(n=40) ... ")
    todo = make_process_list(range(1, 40), fibonacci, return_value=True)
    td1 = time()
    print(execute_process(todo, 10, return_value=True, loop=True))
    print("Time (loop): ", time() - td1)
    td1 = time()
    print(execute_process(todo, 10, return_value=True))
    print("Time (pool): ", time() - td1)


    #
    #
    #
    #
    #   SEARCH PRIMES
    #
    #
    #

    def isprime(n):
        """Returns True if n is prime and False otherwise"""
        if not isinstance(n, int):
            raise TypeError("argument passed to is_prime is not of 'int' type")
        if n < 2:
            return False
        if n == 2:
            return True
        max = int(np.ceil(np.sqrt(n)))
        i = 2
        while i <= max:
            if n % i == 0:
                return False
            i += 1
        return True


    def sum_primes(n):
        """Calculates sum of all primes below given integer n"""
        return sum([x for x in range(2, n) if isprime(x)])


    inputs = range(100000, 1000000, 100000)
    print("Prime Search ... ")
    todo = make_process_list(inputs, sum_primes, return_value=True)
    td1 = time()
    print(execute_process(todo, 10, return_value=True, loop=True))
    print("Time (loop): ", time() - td1)
    td1 = time()
    print(execute_process(todo, 10, return_value=True))
    print("Time (pool): ", time() - td1)
