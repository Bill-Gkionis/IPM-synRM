"""Stand-in programs for the femm_queue tests of the pool workers' jobs. Not a test: each runs in a process of its own, because what it shows would otherwise
happen to pytest itself (a job's limits would hold pytest and every child it starts; an ended owner would be the test run).

    python -I stand_in_pool.py owner   <lock folder> <folder>   own a pool of two workers (no FEMM) whose initializer keeps starting processes, as a worker that
                                                               is opening a FEMM does, and sleep inside it; announce itself first (owner_<pid>.txt)
    python -I stand_in_pool.py parent  <folder>                 the control: no pool, no job. Start one stand-in grandchild below this process and sleep
    python -I stand_in_pool.py no_job  <result file>            run one pool worker's main function here, with femm_queue._Job made to raise (a worker whose job
                                                               cannot be made), and femm_session.open_femm and an initializer that only record that they were
                                                               called; write what the worker told its owner, and what was called, as JSON
    python -I stand_in_pool.py watch   <result file>            run femm_queue._watch_owner() here with an owner that is already gone, and an _end_tree that only
                                                               records its calls; write what the watch thread did (its exit code, and whether it walked a tree)
    python -I stand_in_pool.py join_twice   <result file>       call femm_queue._join_own_job twice here (the first puts this process into a job of its own); write
                                                               what the second call said, and whether the first job was kept
    python -I stand_in_pool.py assign_fails <result file>       call femm_queue._join_own_job here with the job made to refuse this process; write what it said, whether the
                                                               job was closed, and whether this process is in a job
    python -I stand_in_pool.py nested_owner <lock folder> <folder>   announce itself (owner_<pid>.txt), wait for the file <folder>/go (the test puts this process into a job
                                                               meanwhile), then open a pool of two workers (no FEMM; each starts a grandchild) and sleep inside it: the
                                                               workers' own jobs nest under the owner's. Writes <folder>/pool_open.json when the pool is open

-I (isolated mode) keeps the environment's path out, so the script says where femm_queue and queue_test_tools are. A stand-in lives
five minutes at most (the tests end it long before).
"""
import json
import os
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent), str(HERE)]


def owner(lock_dir: str, folder: str) -> None:
    import femm_queue
    import queue_test_tools as qt
    qt.announce(folder, "owner")
    with femm_queue.femm_pool("stand-in owner", workers=2, femm=False, initializer=qt.spawn_grandchildren, initargs=(folder,), lock_dir=Path(lock_dir)):
        time.sleep(300)


def parent(folder: str) -> None:
    import queue_test_tools as qt
    qt.start_grandchild(folder, name="parent")                   # announces this process, then starts the grandchild
    time.sleep(300)


def no_job(result: str) -> None:
    import multiprocessing

    import femm_queue
    import femm_session
    called = {"open_femm": False, "initializer": False}

    def refuse(*arguments, **keywords):
        raise RuntimeError("cannot create a job object for the worker (a stand-in refusal)")

    def open_femm(*arguments, **keywords):
        called["open_femm"] = True

    def initializer():
        called["initializer"] = True
    femm_queue._Job = refuse
    femm_session.open_femm = open_femm
    owner_end, worker_end = multiprocessing.Pipe(duplex=True)
    owner_end.send(None)                                          # "finish up and quit", already waiting: a worker that did start then ends instead of waiting for tasks
    femm_queue._worker_main(worker_end, "below_normal", True, initializer, ())
    message = owner_end.recv()                                    # the first thing the worker said
    Path(result).write_text(json.dumps({"kind": message[0], "text": message[1] if len(message) > 1 else None, **called}), encoding="utf-8")


def watch(result: str) -> None:
    import multiprocessing

    import femm_queue
    calls: list = []
    real_exit = os._exit

    class GoneOwner:                                              # what multiprocessing.parent_process() gives: a process whose join() returns at once
        def join(self):
            return None

    def fake_exit(code):
        Path(result).write_text(json.dumps({"exit_code": code, "end_tree_calls": calls}), encoding="utf-8")
        real_exit(0)
    multiprocessing.parent_process = lambda: GoneOwner()
    femm_queue._end_tree = lambda pid: (calls.append(pid), (0, ""))[1]
    os._exit = fake_exit
    femm_queue._watch_owner()
    time.sleep(30)                                                # the watch thread ends this process long before


def in_any_job() -> bool:
    """Is THIS process in a job (any job)? IsProcessInJob with no job handle asks exactly that."""
    import ctypes
    from ctypes import wintypes
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    k.GetCurrentProcess.restype = wintypes.HANDLE
    k.IsProcessInJob.argtypes = [wintypes.HANDLE, wintypes.HANDLE, ctypes.POINTER(wintypes.BOOL)]
    k.IsProcessInJob.restype = wintypes.BOOL
    inside = wintypes.BOOL()
    if not k.IsProcessInJob(k.GetCurrentProcess(), None, ctypes.byref(inside)):
        raise ctypes.WinError(ctypes.get_last_error())
    return bool(inside.value)


def join_twice(result: str) -> None:
    import femm_queue
    before = in_any_job()
    femm_queue._join_own_job("below_normal")                      # the first call: this process is in a job of its own now
    first = femm_queue._OWN_JOB
    after_first = in_any_job()
    try:
        femm_queue._join_own_job("below_normal")
        second = None
    except RuntimeError as error:
        second = str(error)
    Path(result).write_text(json.dumps({"in_a_job_before": before, "in_a_job_after_the_first": after_first, "second": second,
                                        "first_job_kept": femm_queue._OWN_JOB is first, "first_job_handle_open": bool(first.handle)}), encoding="utf-8")


def assign_fails(result: str) -> None:
    import femm_queue
    closed = []
    real_close = femm_queue._Job.close

    def refuse(self, pid):
        raise RuntimeError(f"cannot put {self.who} into its job: [WinError 5] Access is denied. (a stand-in refusal)")

    def close(self):
        closed.append(self.who)
        real_close(self)
    femm_queue._Job.assign, femm_queue._Job.close = refuse, close
    before = in_any_job()
    try:
        femm_queue._join_own_job("below_normal")
        said = None
    except RuntimeError as error:
        said = str(error)
    Path(result).write_text(json.dumps({"said": said, "closed": closed, "own_job_kept": femm_queue._OWN_JOB is not None, "in_a_job_before": before,
                                        "in_a_job_after": in_any_job()}), encoding="utf-8")


def nested_owner(lock_dir: str, folder: str) -> None:
    import femm_queue
    import queue_test_tools as qt
    qt.announce(folder, "owner")
    deadline = time.monotonic() + 120
    while not Path(folder, "go").exists() and time.monotonic() < deadline:
        time.sleep(0.02)
    with femm_queue.femm_pool("nested owner", workers=2, femm=False, initializer=qt.start_grandchild, initargs=(folder,), lock_dir=Path(lock_dir)) as pool:
        facts = pool.map(qt.own_job_facts, [0, 1, 2, 3])          # (worker pid, has a job of its own, is in it, its handle is inheritable)
        Path(folder, "pool_open.json").write_text(json.dumps({"workers": pool.pids, "facts": facts}), encoding="utf-8")
        time.sleep(300)


if __name__ == "__main__":
    {"owner": owner, "parent": parent, "no_job": no_job, "watch": watch, "join_twice": join_twice, "assign_fails": assign_fails,
     "nested_owner": nested_owner}[sys.argv[1]](*sys.argv[2:])
