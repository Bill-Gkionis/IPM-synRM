"""Tests for femm_queue: the machine-wide solver slot, the FEMM pool, and run_solver.

None of these starts a real solver, and only two start FEMM (both are marked `femm`, so `-m "not femm"` leaves them out). The "other
run" is a stand-in process that only holds the slot, and the pool's tasks are small functions that return, hang, crash or start a
child (tests/queue_test_tools.py). Every test works in a lock folder of its own (the `lock_dir` fixture), never in the real
out/.femm: a real run may be holding that slot while the tests run. The pools are told femm=False, which starts no FEMM in their
workers; the live check with real FEMMs, two pools and a solver at once, is check_queue.py.

Run from P:   uv run pytest tests/test_femm_queue.py -q       (about 40 s on 2026-11-06 without the two `femm` tests: it starts a few dozen small Python processes)

The seven tests the plan names are here under its names: test_second_caller_waits_then_runs, test_busy_names_the_holder,
test_killed_holder_frees_the_queue, test_holder_is_none_when_free, test_stuck_task_is_killed_and_retried,
test_ctrl_c_while_waiting_leaves_no_lock and test_run_solver_kills_on_timeout. The rest pin what the module promises beyond them.
"""
import json
import os
import re
import subprocess
import sys
import threading
import time
import types
from datetime import datetime
from pathlib import Path

import pytest

import femm_queue
import femm_test_tools as tools
import paths
import queue_test_tools as qt

P = Path(__file__).resolve().parent.parent
BELOW_NORMAL = 0x4000                      # BELOW_NORMAL_PRIORITY_CLASS
IDLE = 0x40                                # IDLE_PRIORITY_CLASS


@pytest.fixture
def lock_dir(tmp_path):
    """The lock folder of this test alone. No test touches the real slot in out/.femm."""
    return tmp_path / "slot"


def write_config(folder: Path, workers: int = 8, priority: str = "below_normal") -> Path:
    """A config.toml of its own, so a test can change the limits (paths.CONFIG is read at call time) and the lock folder."""
    config = folder / "config.toml"
    config.write_text(f"[paths]\nout = '{folder / 'out'}'\n[femm]\nworkers = {workers}\npriority = '{priority}'\n", encoding="utf-8")
    return config


# ======================================================================================= the slot: who holds it, and who waits
def test_holder_is_none_when_free(lock_dir):
    """Nobody holds the slot: holder() says None. Before the first use, after a use, and with a dead run's note left behind.

    The last case is the one that matters. holder.json is a note for people; the lock is the truth. A run that was killed
    leaves its note behind, and holder() must not report a run that no longer exists.
    """
    assert femm_queue.holder(lock_dir) is None
    assert not lock_dir.exists(), "asking who holds the slot must not create anything"
    with femm_queue.femm_slot("brief", lock_dir=lock_dir):
        assert femm_queue.holder(lock_dir)["label"] == "brief"
    assert femm_queue.holder(lock_dir) is None
    (lock_dir / "holder.json").write_text(json.dumps({"label": "ghost", "pid": 4_000_000_000, "since": "2020-01-01T00:00:00"}))
    assert femm_queue.holder(lock_dir) is None


def test_a_free_slot_is_taken_at_once_and_everything_is_given_back(lock_dir):
    """A free slot is taken without waiting (wait=False would raise otherwise); the note is complete; only the lock file stays."""
    with femm_queue.solver_slot("bench B4 ccx", wait=False, lock_dir=lock_dir) as info:
        note = json.loads((lock_dir / "holder.json").read_text())
        assert note == info and set(info) == {"label", "pid", "since"}
        assert info["label"] == "bench B4 ccx" and info["pid"] == os.getpid()
        datetime.fromisoformat(info["since"])              # a time the waiting line can print as HH:MM
    assert sorted(p.name for p in lock_dir.iterdir()) == ["solver.lock"], "a ticket or a note was left behind"


def test_femm_slot_is_the_solver_slot():
    """R10: one lock for every heavy run. femm_slot is the first name of the same thing, not a second lock."""
    assert femm_queue.femm_slot is femm_queue.solver_slot


def test_busy_names_the_holder(lock_dir):
    """wait=False on a taken slot raises QueueBusy, and the message says who: the label, the pid and the time of day."""
    assert issubclass(femm_queue.QueueBusy, RuntimeError)
    with qt.Holder(lock_dir, "bench B4 ccx") as first:
        with pytest.raises(femm_queue.QueueBusy) as busy:
            with femm_queue.femm_slot("me", wait=False, lock_dir=lock_dir):
                pytest.fail("took a slot that another process holds")
        message = str(busy.value)
        assert "'bench B4 ccx'" in message and f"(pid {first.pid})" in message
        assert re.search(r"since \d\d:\d\d\b", message), message
        info = femm_queue.holder(lock_dir)
        assert (info["label"], info["pid"]) == ("bench B4 ccx", first.pid)
        assert busy.value.holder == info


def test_second_caller_waits_then_runs(lock_dir, tmp_path, capsys):
    """The second caller prints a waiting line every poll, does not run while the first holds, and runs once it lets go.

    "Does not overlap" is shown from the clock: the first holder writes the time it is about to give the slot back, the
    second one the time it got it, and the second is never earlier. The waiter is a thread here: a lock is per open file, so a
    second handle in this very process competes with the other process's exactly as a second terminal would.
    """
    log = tmp_path / "first.log"
    ran = {}

    def second():
        with femm_queue.femm_slot("second", poll_s=0.2, lock_dir=lock_dir):
            ran["at"] = time.time()
            ran["holder"] = femm_queue.holder(lock_dir)

    with qt.Holder(lock_dir, "first", log=log) as first:
        waiter = threading.Thread(target=second, daemon=True)
        waiter.start()
        time.sleep(1.2)                                        # six polls of 0.2 s
        assert waiter.is_alive() and "at" not in ran, "the second caller ran while the first held the slot"
        waiting = [line for line in capsys.readouterr().out.splitlines() if line.startswith("waiting for the solver slot")]
        assert len(waiting) >= 2, "one line per poll was expected"
        assert re.fullmatch(rf"waiting for the solver slot: held by 'first' \(pid {first.pid}\) since \d\d:\d\d", waiting[0]), waiting[0]
        first.release()
        waiter.join(30)
        assert not waiter.is_alive(), "the second caller never got the slot after the first let go"
    assert ran["holder"]["label"] == "second" and ran["holder"]["pid"] == os.getpid()
    gave_back = float(log.read_text().splitlines()[-1].split()[1])
    assert ran["at"] >= gave_back, "the two runs overlapped"
    assert femm_queue.holder(lock_dir) is None


def test_killed_holder_frees_the_queue(lock_dir):
    """A holder that is killed (no `finally`, no clean-up, nothing) frees the slot: the operating system lets go of its lock.

    Its note is still there afterwards, and holder() ignores it. That is the whole argument for a lock over a PID file.
    """
    with qt.Holder(lock_dir, "doomed") as doomed:
        assert femm_queue.holder(lock_dir)["pid"] == doomed.pid
        with pytest.raises(femm_queue.QueueBusy):
            with femm_queue.femm_slot("early", wait=False, lock_dir=lock_dir):
                pass
        doomed.kill()
    deadline = time.monotonic() + 5                            # the OS needs no more than an instant; this only absorbs a slow PC
    while femm_queue.holder(lock_dir) is not None and time.monotonic() < deadline:
        time.sleep(0.05)
    assert (lock_dir / "holder.json").is_file(), "the killed run cannot have cleaned up"
    assert femm_queue.holder(lock_dir) is None
    with femm_queue.femm_slot("next", wait=False, lock_dir=lock_dir):
        assert femm_queue.holder(lock_dir)["label"] == "next"


def test_ctrl_c_while_waiting_leaves_no_lock(lock_dir):
    """A caller that is waiting and gets a Ctrl-C ends with KeyboardInterrupt and leaves nothing behind.

    The check that matters is the holder's note: a clean-up that deleted holder.json unconditionally would erase the note of
    the run that really holds the slot, and the next waiting caller would print nothing useful. Afterwards, once the real
    holder lets go, the slot is free at once and the folder holds only the lock file.
    """
    outcome = qt.CTX.Queue()
    with qt.Holder(lock_dir, "the real holder") as real:
        waiter = qt.CTX.Process(target=qt.wait_then_get_interrupted, args=(str(lock_dir), "impatient", 0.8, outcome))
        waiter.start()
        assert outcome.get(timeout=60) == "interrupted"
        waiter.join(30)
        info = femm_queue.holder(lock_dir)
        assert (info["label"], info["pid"]) == ("the real holder", real.pid), "the interrupted waiter disturbed the holder's note"
        real.release()
    assert femm_queue.holder(lock_dir) is None
    with femm_queue.femm_slot("after", wait=False, lock_dir=lock_dir):
        pass
    assert sorted(p.name for p in lock_dir.iterdir()) == ["solver.lock"]


def test_an_error_inside_the_slot_still_gives_it_back(lock_dir):
    """Whatever ends the `with` (an exception, a Ctrl-C), the slot is given back and the note removed."""
    for failure in (RuntimeError("the solver crashed"), KeyboardInterrupt()):
        with pytest.raises(type(failure)):
            with femm_queue.femm_slot("crashing", wait=False, lock_dir=lock_dir):
                raise failure
        assert femm_queue.holder(lock_dir) is None
        assert sorted(p.name for p in lock_dir.iterdir()) == ["solver.lock"]


def test_a_second_slot_in_the_same_thread_is_refused_not_waited_for(lock_dir):
    """A run that asks for the slot it already holds would wait for itself for ever. It gets an error that says so instead.

    wait=False below only so that a missing guard fails this test (QueueBusy, naming this very process as the holder) and does
    not hang it: the guard sits in front of the waiting, so wait=True takes the same road.
    """
    with femm_queue.femm_slot("outer", lock_dir=lock_dir):
        with pytest.raises(RuntimeError, match="already holds the solver slot") as refused:
            with femm_queue.femm_slot("inner", wait=False, lock_dir=lock_dir):
                pass
        assert not isinstance(refused.value, femm_queue.QueueBusy), "nested use is a mistake, not a busy slot"
        assert femm_queue.holder(lock_dir)["label"] == "outer"
    assert femm_queue.holder(lock_dir) is None


def test_the_lock_folder_and_the_limits_come_from_the_config(tmp_path, monkeypatch):
    """Default lock folder = [paths] out + ".femm"; the worker and thread limit = [femm] workers; the priority = [femm] priority.

    The refusals take no slot and start nothing; 'normal' priority is refused, because rule R10 is "below normal".
    """
    config = write_config(tmp_path, workers=2)
    monkeypatch.setattr(paths, "CONFIG", config)
    assert femm_queue.holder() is None
    with femm_queue.solver_slot("from the config", wait=False):
        assert (tmp_path / "out" / ".femm" / "holder.json").is_file()
        assert femm_queue.holder()["label"] == "from the config"
    assert femm_queue.holder() is None
    with pytest.raises(ValueError, match=r"workers=3.*between 1 and 2"):
        with femm_queue.femm_pool("too many", workers=3, femm=False):
            pass
    quick = [sys.executable, "-c", "pass"]
    with pytest.raises(ValueError, match=r"threads=3.*between 1 and 2"):
        femm_queue.run_solver(quick, label="too many", cwd=tmp_path, timeout_s=5, threads=3)
    write_config(tmp_path, priority="normal")
    with pytest.raises(ValueError, match="priority"):
        femm_queue.run_solver(quick, label="normal priority", cwd=tmp_path, timeout_s=5)
    with pytest.raises(ValueError, match="priority"):
        with femm_queue.femm_pool("normal priority", workers=1, femm=False):
            pass
    assert femm_queue.holder() is None


def test_idle_priority_is_honoured_when_the_config_asks_for_it(tmp_path, lock_dir, monkeypatch):
    """[femm] priority = 'idle' is the one other class the module accepts (lower still than below normal): a solver and a worker get it."""
    monkeypatch.setattr(paths, "CONFIG", write_config(tmp_path, priority="idle"))
    code = ("import ctypes\n"
            "k = ctypes.WinDLL('kernel32')\n"
            "k.GetCurrentProcess.restype = ctypes.c_void_p\n"
            "k.GetPriorityClass.argtypes = [ctypes.c_void_p]\n"
            "print(k.GetPriorityClass(k.GetCurrentProcess()))\n")
    run = femm_queue.run_solver([sys.executable, "-c", code], label="idle solver", cwd=tmp_path, timeout_s=60, lock_dir=lock_dir)
    assert run.stdout.split() == [str(IDLE)], run
    with femm_queue.femm_pool("idle pool", workers=1, femm=False, lock_dir=lock_dir) as pool:
        assert [tools.priority_class(pid) for pid in pool.pids] == [IDLE]


# ======================================================================================= the pool
def test_the_pool_runs_tasks_on_its_own_workers_at_below_normal_priority(tmp_path, lock_dir, monkeypatch):
    """Results come back in the order of the items, from at most `workers` processes of the pool, below normal priority.

    workers=None takes the configured number (2 here). The pool holds the slot while it is open, a task cannot take that slot
    again (it is refused, not left waiting), and when the `with` ends every worker is gone and the slot is free.
    """
    monkeypatch.setattr(paths, "CONFIG", write_config(tmp_path, workers=2))
    with femm_queue.femm_pool("pool basics", femm=False, lock_dir=lock_dir) as pool:
        pids = list(pool.pids)
        assert pool.workers == 2 and len(pids) == 2 and os.getpid() not in pids
        results = pool.map(qt.echo_pid, list(range(7)))
        assert [item for item, _ in results] == list(range(7))
        assert {pid for _, pid in results} <= set(pids)
        assert [tools.priority_class(pid) for pid in pids] == [BELOW_NORMAL, BELOW_NORMAL]
        assert pool.map(qt.try_to_take_the_slot, [str(lock_dir)])[0].startswith("REFUSED"), "a worker may not take its pool's slot"
        assert femm_queue.holder(lock_dir)["label"] == "pool basics"
        assert pool.map(qt.echo_pid, []) == []
    assert qt.wait_until_gone(pids) == [], "a worker was still running after the pool closed"
    assert femm_queue.holder(lock_dir) is None


def test_stuck_task_is_killed_and_retried(tmp_path, lock_dir):
    """A task over task_timeout_s has its worker's whole process tree ended; the task is run once more by another worker.

    Item 1 hangs on its first attempt only (and starts a child below the worker, like FEMM's fkn.exe). The worker that
    hung is gone, its child is gone, a different process made the second attempt, the pool is back at full strength, and
    every item has its result, in order.
    """
    items = [(str(tmp_path), index, "first" if index == 1 else "never") for index in range(4)]
    started = time.monotonic()
    with femm_queue.femm_pool("stuck once", workers=2, task_timeout_s=3, femm=False, lock_dir=lock_dir) as pool:
        before = set(pool.pids)
        results = pool.map(qt.maybe_stuck, items)
        after = set(pool.pids)
    elapsed = time.monotonic() - started
    assert [index for index, _ in results] == [0, 1, 2, 3]
    first_try, second_try = [int(pid) for pid in (tmp_path / "item1.pids").read_text().split()]
    assert first_try != second_try, "the retry ran in the worker that was stuck"
    assert first_try in before and first_try not in after and len(after) == 2, "the stuck worker was not replaced"
    child = int((tmp_path / "item1.children").read_text())
    assert qt.wait_until_gone([first_try, child]) == [], "the stuck worker or the child below it was left running"
    assert 3 <= elapsed < 60


def test_a_task_stuck_twice_is_a_timeout_and_only_the_pools_own_processes_are_ended(tmp_path, lock_dir):
    """Stuck on the retry too: the item's result is {"status": "timeout"}, and the run goes on to the other items.

    A bystander (a python.exe that is not the pool's, started first) must survive both kills: the pool ends processes by
    the pid it started, never by name.
    """
    bystander = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
    try:
        items = [(str(tmp_path), index, "always" if index == 1 else "never") for index in range(3)]
        with femm_queue.femm_pool("stuck twice", workers=2, task_timeout_s=2.5, femm=False, lock_dir=lock_dir) as pool:
            results = pool.map(qt.maybe_stuck, items)
        assert results[1] == {"status": "timeout"}
        assert [r[0] for r in (results[0], results[2])] == [0, 2], "the other items were lost"
        pids = [int(p) for p in (tmp_path / "item1.pids").read_text().split()]
        children = [int(p) for p in (tmp_path / "item1.children").read_text().split()]
        assert len(pids) == 2 and pids[0] != pids[1], "the item was tried once, or twice in the same worker"
        assert qt.wait_until_gone(pids + children) == [], "a stuck worker or its child was left running"
        assert bystander.poll() is None, "a process that the pool did not start was ended"
        assert femm_queue.holder(lock_dir) is None
    finally:
        bystander.kill()


def test_a_task_that_raises_is_raised_in_the_caller_and_the_slot_is_given_back(lock_dir):
    """An exception in a task is not a timeout: it comes out of map() as the same exception, with the worker's traceback."""
    with pytest.raises(ValueError, match="boom 3") as raised:
        with femm_queue.femm_pool("raising", workers=1, femm=False, lock_dir=lock_dir) as pool:
            pids = list(pool.pids)
            pool.map(qt.fail_on_three, [1, 2, 3, 4])
    assert "fail_on_three" in str(raised.value.__cause__), "the worker's traceback should travel with the exception"
    assert qt.wait_until_gone(pids) == []
    assert femm_queue.holder(lock_dir) is None


def test_a_task_that_raises_an_oserror_is_an_error_not_a_dead_worker(tmp_path, lock_dir):
    """A FileNotFoundError in a task (a model file that is not there) is the task's error, raised in the caller, at the first
    attempt. It must not be taken for a worker that died (a broken pipe is an OSError too) and be retried into {"status": "died"}."""
    with pytest.raises(FileNotFoundError, match="no model file motor.fem") as raised:
        with femm_queue.femm_pool("missing file", workers=1, femm=False, lock_dir=lock_dir) as pool:
            pool.map(qt.missing_file, [(str(tmp_path), "motor.fem")])
    assert "missing_file" in str(raised.value.__cause__)
    assert len((tmp_path / "attempts.txt").read_text().split()) == 1, "the task was tried again"
    assert femm_queue.holder(lock_dir) is None


def test_a_worker_that_dies_is_replaced_and_its_task_tried_once_more(tmp_path, lock_dir):
    """A worker that ends by itself (here: os._exit, a crash) is handled like a timeout: replaced, its task tried once more,
    and after a second crash the item's result is {"status": "died"}. The pool is back at full strength afterwards."""
    items = [(str(tmp_path), 0, "never"), (str(tmp_path), 1, "first"), (str(tmp_path), 2, "always")]
    with femm_queue.femm_pool("dying", workers=2, task_timeout_s=60, femm=False, lock_dir=lock_dir) as pool:
        results = pool.map(qt.maybe_die, items)
        assert len(pool.pids) == 2
    assert results[0][0] == 0 and results[1][0] == 1 and results[2] == {"status": "died"}
    for index in (1, 2):
        pids = (tmp_path / f"item{index}.pids").read_text().split()
        assert len(pids) == 2 and pids[0] != pids[1], f"item {index} was not tried once on a fresh worker"
    assert femm_queue.holder(lock_dir) is None


def test_a_function_a_worker_cannot_import_is_an_error_not_a_dead_worker(lock_dir):
    """The classic spawn mistake: a function that the owner can pickle but a new Python cannot find.

    The module below exists only in this process (sys.modules), so a worker fails to unpickle the function. That has to come out
    of map() as the worker's own ModuleNotFoundError with its traceback, and not as a worker that "died twice". A lambda is
    refused in the caller, before any worker is bothered.
    """
    module = types.ModuleType("only_in_the_owner")
    exec("def double(x):\n    return 2 * x\n", module.__dict__)
    sys.modules[module.__name__] = module
    try:
        with pytest.raises(ModuleNotFoundError, match="only_in_the_owner") as raised:
            with femm_queue.femm_pool("unimportable", workers=1, femm=False, lock_dir=lock_dir) as pool:
                pool.map(module.double, [1, 2])
        assert "only_in_the_owner" in str(raised.value.__cause__), "the worker's traceback should travel with the error"
        with femm_queue.femm_pool("lambda", workers=1, femm=False, lock_dir=lock_dir) as pool:
            with pytest.raises(TypeError, match="top level"):
                pool.map(lambda item: item, [1])
    finally:
        del sys.modules[module.__name__]
    assert femm_queue.holder(lock_dir) is None


def test_a_worker_that_hangs_or_fails_while_starting_is_reported_not_waited_for(tmp_path, lock_dir):
    """The start-up is under a timeout too: a FEMM that never answers cannot hold the slot for ever.

    This is the net femm_session's docstring points to: open_femm() has no timeout of its own, so the pool's start timeout ends a
    worker that is stuck in it. An initializer that fails is reported with its traceback.
    """
    pid_file = tmp_path / "hung.pid"
    started = time.monotonic()
    with pytest.raises(RuntimeError, match="did not start"):
        with femm_queue.femm_pool("hang at start", workers=1, femm=False, start_timeout_s=2, initializer=qt.hang_forever,
                                  initargs=(str(pid_file),), lock_dir=lock_dir):
            pytest.fail("the pool opened although its worker never started")
    assert time.monotonic() - started < 60
    assert qt.wait_until_gone([int(pid_file.read_text())]) == [], "the worker that hung while starting was left running"
    assert femm_queue.holder(lock_dir) is None
    with pytest.raises(RuntimeError, match="failed while starting") as failed:
        with femm_queue.femm_pool("fail at start", workers=1, femm=False, initializer=qt.fail_to_start, lock_dir=lock_dir):
            pytest.fail("the pool opened although its worker's initializer failed")
    assert "no FEMM here" in str(failed.value)
    assert femm_queue.holder(lock_dir) is None
    with pytest.raises(RuntimeError, match=r"died while starting \(exit code 5\)"):
        with femm_queue.femm_pool("crash at start", workers=1, femm=False, initializer=qt.die_while_starting, lock_dir=lock_dir):
            pytest.fail("the pool opened although its worker died while starting")
    assert femm_queue.holder(lock_dir) is None


def test_workers_end_when_the_process_that_owns_the_pool_is_killed(lock_dir):
    """Kill the pool's owner (no clean-up at all): its workers end too, so nothing goes on computing after the slot is free.

    The workers are in the middle of a task that never ends (an idle worker would end by itself, as its pipe closes). What ends
    a busy one is the thread every worker starts to watch its owner. On Windows it ends the worker (os._exit) and the worker's own
    job ends everything below it (its FEMM); elsewhere it ends the worker's process group, which takes the same with it.
    """
    pids_out = qt.CTX.Queue()
    owner = qt.CTX.Process(target=qt.pool_owner, args=(str(lock_dir), pids_out))
    owner.start()
    workers = []                                                # handles, opened as soon as the pool reports its workers
    try:
        reported = pids_out.get(timeout=60)                     # (pid, creation time) of each worker, taken where the pool holds them
        workers = [qt.Owned(pid, created) for pid, created in reported]     # running, and the pool's: the handles are ours
        assert len(workers) == 2
        owner.kill()
        owner.join(10)
        left = [w.pid for w in workers if not w.wait_gone(20)]
        assert left == [], f"workers {left} outlived the process that owned their pool"
        assert femm_queue.holder(lock_dir) is None
    finally:
        if owner.is_alive():
            owner.kill()
        for w in workers:                                       # a failed test leaves nothing: ended through the handle, never by pid
            w.end()
            w.close()


# ======================================================================================= the pool workers' jobs (the finding of a timeout trial of B9)
# A timeout trial of B9 on a real bench tree (2026-11-06) ended a pool's owner 0.1 s into the pool's life. The owner's eight workers went, but one FEMM did not: it was created by its
# worker's main thread (femm_session.open_femm) a moment after the worker's watch thread had walked the worker's tree with taskkill /T, and before the worker was gone, so
# the walk did not have it. The fix has no walk to be early for: on Windows every pool worker joins a kill-on-close job of its own before it starts any process, and the
# system ends whatever is in it when the worker's last handle closes, however the worker ends. These tests show it with stand-in processes (a grandchild that sleeps two
# minutes), never with the test process itself in a job (its limits would then hold pytest and every child it starts): a worker, or an owner, or a plain parent is
# ended through a handle the test opens on it (qt.Owned: pid and creation time, checked), and whatever a test leaves is ended through its own handle, by pid and
# creation time, in a finally.
STAND_IN_POOL = str(Path(__file__).with_name("stand_in_pool.py"))      # the stand-in programs of these tests (a script of its own)


def pid_and_time(path) -> tuple[int, int]:
    """(pid, creation time) from a file that a stand-in announced itself in (qt.announce)."""
    pid, created = (int(part) for part in Path(path).read_text().split())
    return pid, created


def grandchildren_written_by(folder, worker_pid=None) -> list[tuple[int, int]]:
    """Every (pid, creation time) that workers wrote in <folder>/grandchildren_<worker pid>.txt as they started a grandchild (one worker's, or all)."""
    found = []
    for path in sorted(Path(folder).glob("grandchildren_*.txt" if worker_pid is None else f"grandchildren_{worker_pid}.txt")):
        for line in path.read_text().splitlines():
            parts = line.split()
            if len(parts) == 2:
                found.append((int(parts[0]), int(parts[1])))
    return found


def end_each_written_down(processes) -> None:
    """The tests' clean-up: each (pid, creation time) that is still running is ended through one handle opened on it, and only if the time read through the handle
    is the one that was written down (a number given to a stranger since is refused). So a test that fails leaves nothing of its own running, and ends nothing else."""
    for pid, created in processes:
        try:
            owned = qt.Owned(pid, created)
        except ProcessLookupError:
            continue
        owned.end()
        owned.close()


def test_a_worker_that_is_ended_with_no_clean_up_takes_what_it_started_with_it(tmp_path, lock_dir):
    """The worker's death takes its children. Two workers, each of which starts a stand-in grandchild in its initializer (before it says it is ready). One worker is ended
    with TerminateProcess through a handle: no tree walk, no watch thread, no taskkill, nothing but that call. Its grandchild must be gone within 5 s, and the other
    worker's grandchild must still run (each worker has a job of its own: one job for the pool would end both).

    Before the change the grandchild survives (nothing walks the tree of a worker that is ended from outside), and the finally ends it, so nothing is left either way."""
    handles = []
    try:
        with femm_queue.femm_pool("worker death", workers=2, femm=False, initializer=qt.start_grandchild, initargs=(str(tmp_path),), lock_dir=lock_dir) as pool:
            workers = {pid: qt.creation_time_of(pid) for pid in pool.pids}                  # taken here, where the pool holds the workers: each number is certainly that worker's
            assert len(workers) == 2
            first, second = list(workers)
            grandchildren = {}
            for pid in workers:
                (g_pid, g_created), = grandchildren_written_by(tmp_path, pid)
                grandchildren[pid] = qt.Owned(g_pid, g_created)
                handles.append(grandchildren[pid])
            victim = qt.Owned(first, workers[first])
            handles.append(victim)
            assert grandchildren[first].running() and grandchildren[second].running()
            victim.end()                                                                      # TerminateProcess through the handle: the worker ends with no clean-up at all
            assert victim.wait_gone(10), "the worker was not ended"
            assert grandchildren[first].wait_gone(5), "the grandchild of the ended worker is still running 5 s after its worker ended"
            assert grandchildren[second].running(), "the other worker's grandchild went too: the jobs are one for each worker"
    finally:
        for handle in handles:
            handle.end()                                                                      # (a handle on a process that has ended does nothing)
            handle.close()
        end_each_written_down(grandchildren_written_by(tmp_path))                             # one that was written down but got no handle (the test failed before it): ended too


def test_a_worker_that_hangs_past_the_timeout_is_ended_through_the_pools_handle_and_its_job_ends_the_rest_with_no_taskkill(tmp_path, lock_dir, monkeypatch):
    """Final review, Important 2. A task over task_timeout_s has its worker ended by the pool. On Windows that is TerminateProcess through the handle the pool holds on the worker
    (Process.kill), and the worker's own job ends what the worker started. No taskkill /T is started: it walks the tree through each process's RECORDED parent number, so it could
    take an orphan whose dead parent's number the worker now holds, a process this run did not start.

    Two workers, each of which starts a stand-in grandchild in its initializer. Item 1 hangs on its first attempt only. The owner (this very process) is watched: it starts no program
    with subprocess.run and calls _end_tree never. The hung worker goes, its grandchild is gone within 5 s of it, the other worker's grandchild still runs (a job for each worker, not
    one for the pool), the pool is back at two workers (a fresh one took the hung worker's place), and item 1 was run again, by another worker, and returned.

    Before the change a taskkill was started (the grandchild went too: the walk found it), so what turns this test red is the taskkill. Every stand-in is ended in the finally, by
    pid and creation time, through a handle."""
    started, walked = [], []
    real_run, real_end_tree = subprocess.run, femm_queue._end_tree
    monkeypatch.setattr(subprocess, "run", lambda *args, **kwargs: (started.append(args[0] if args else kwargs.get("args")), real_run(*args, **kwargs))[1])
    monkeypatch.setattr(femm_queue, "_end_tree", lambda pid: (walked.append(pid), real_end_tree(pid))[1])
    handles, outcome = [], {}
    try:
        with femm_queue.femm_pool("hung worker", workers=2, femm=False, task_timeout_s=3, initializer=qt.start_grandchild, initargs=(str(tmp_path),),
                                  lock_dir=lock_dir) as pool:
            items = [(str(tmp_path), 0, "never"), (str(tmp_path), 1, "first")]
            runner = threading.Thread(target=lambda: outcome.update(results=pool.map(qt.hang_the_first_time, items)), daemon=True)
            try:
                workers = {pid: qt.creation_time_of(pid) for pid in pool.pids}              # taken where the pool holds the workers: each number is certainly that worker's
                assert len(workers) == 2
                owned, grandchild = {}, {}
                for pid, created in workers.items():
                    owned[pid] = qt.Owned(pid, created)
                    (g_pid, g_created), = grandchildren_written_by(tmp_path, pid)
                    grandchild[pid] = qt.Owned(g_pid, g_created)
                    handles += [owned[pid], grandchild[pid]]
                runner.start()
                deadline = time.monotonic() + 60
                while not (tmp_path / "item1.pids").exists() and time.monotonic() < deadline:
                    time.sleep(0.01)
                victim = int((tmp_path / "item1.pids").read_text().split()[0])              # the worker that took item 1 and hangs in it
                assert victim in owned, "item 1 was taken by a process that is not one of the pool's workers"
                survivor = next(pid for pid in owned if pid != victim)
                assert owned[victim].wait_gone(30), "the hung worker was still running 30 s after its task started (the timeout is 3 s)"
                assert grandchild[victim].wait_gone(5), "the grandchild of the ended worker is still running 5 s after its worker went"
                assert started == [] and walked == [], f"the owner started {started} and walked {walked}: a worker is ended through its handle, with no taskkill"
                assert grandchild[survivor].running(), "the other worker's grandchild went too: each worker has a job of its own"
            finally:
                if runner.is_alive():
                    runner.join(90)                                                          # the batch ends by itself within seconds; the pool must not close under it
            assert not runner.is_alive(), "the batch did not finish"
            assert [result[0] for result in outcome["results"]] == [0, 1]
            tries = [int(pid) for pid in (tmp_path / "item1.pids").read_text().split()]
            assert len(tries) == 2 and tries[0] == victim != tries[1], f"item 1 was tried by workers {tries}: the retry has to run in another worker"
            assert outcome["results"][1][1] == tries[1] and len(pool.pids) == 2 and victim not in pool.pids, "the pool did not go on with a fresh worker"
    finally:
        for handle in handles:
            handle.end()                                                                      # (a handle on a process that has ended does nothing)
            handle.close()
        end_each_written_down(grandchildren_written_by(tmp_path))
    assert started == [] and walked == []


def test_the_owner_ended_while_its_workers_are_starting_processes_leaves_nothing_running(tmp_path, lock_dir):
    """The path that trial took. A stand-in owner (a `python -I` script) opens a pool of two workers whose initializers keep starting grandchildren with their main thread, for about
    a second, as a worker does while it opens its FEMM. As soon as both are at it the owner is ended with TerminateProcess through a handle. Every worker, and every
    grandchild either wrote down (they write each one as it starts), must be gone within 5 s.

    Before the change this may pass by timing: the watch thread's taskkill finds what exists when it walks, and what a worker's main thread creates after the walk and before
    the worker is gone survives. Which ones that is, if any, is a matter of milliseconds."""
    launcher = subprocess.Popen([sys.executable, "-I", STAND_IN_POOL, "owner", str(lock_dir), str(tmp_path)], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL)
    handles = []
    try:
        def both_are_at_it() -> bool:
            return bool(list(tmp_path.glob("owner_*.txt"))) and len(list(tmp_path.glob("worker_*.txt"))) == 2 and len(list(tmp_path.glob("grandchildren_*.txt"))) == 2
        deadline = time.monotonic() + 90
        while not both_are_at_it() and time.monotonic() < deadline:
            time.sleep(0.005)
        assert both_are_at_it(), "the stand-in owner's workers did not start in 90 s"
        owner = qt.Owned(*pid_and_time(next(tmp_path.glob("owner_*.txt"))))                   # the interpreter of the stand-in (the launcher in front of it is not the owner)
        workers = [qt.Owned(*pid_and_time(path)) for path in sorted(tmp_path.glob("worker_*.txt"))]
        handles += [owner, *workers]
        owner.end()                                                                           # TerminateProcess: the workers are in the middle of their initializers
        assert owner.wait_gone(10), "the owner was not ended"
        left = [worker.pid for worker in workers if not worker.wait_gone(5)]
        assert left == [], f"workers {left} outlived the process that owned their pool"
        survivors = []
        for pid, created in grandchildren_written_by(tmp_path):                               # the workers are gone: nothing is being added any more
            try:
                grandchild = qt.Owned(pid, created)
            except ProcessLookupError:                                                        # gone already
                continue
            handles.append(grandchild)
            if not grandchild.wait_gone(5):
                survivors.append(pid)
        assert survivors == [], f"{len(survivors)} process(es) that a worker had started were still running 5 s after its owner ended: {survivors[:6]}"
    finally:
        launcher.kill()
        launcher.wait(30)
        for handle in handles:
            handle.end()
            handle.close()
        end_each_written_down(grandchildren_written_by(tmp_path))


def test_without_a_pool_and_its_job_a_process_ended_with_no_clean_up_leaves_what_it_started_running(tmp_path):
    """The control for the two tests above. The same grandchild below a plain stand-in parent: no pool, no job. The parent is ended with TerminateProcess through a handle, and
    the grandchild keeps running (5 s later it still does). So those tests can tell a job from no job; without it they would show nothing. The grandchild is ended by the test."""
    launcher = subprocess.Popen([sys.executable, "-I", STAND_IN_POOL, "parent", str(tmp_path)], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    handles = []
    try:
        deadline = time.monotonic() + 60
        while not (list(tmp_path.glob("parent_*.txt")) and list(tmp_path.glob("grandchildren_*.txt"))) and time.monotonic() < deadline:
            time.sleep(0.01)
        parent = qt.Owned(*pid_and_time(next(tmp_path.glob("parent_*.txt"))))
        (g_pid, g_created), = grandchildren_written_by(tmp_path)
        grandchild = qt.Owned(g_pid, g_created)
        handles += [parent, grandchild]
        parent.end()
        assert parent.wait_gone(10), "the parent was not ended"
        assert not grandchild.wait_gone(5), "the grandchild went with its parent although nothing made it: the control shows no difference, so the tests above show nothing"
    finally:
        launcher.kill()
        launcher.wait(30)
        for handle in handles:
            handle.end()
            handle.close()
        end_each_written_down(grandchildren_written_by(tmp_path))


def test_a_worker_whose_job_cannot_be_made_says_so_and_never_starts_a_femm(tmp_path):
    """The failure path. A worker that cannot join its job does not run unprotected, silently: it tells its owner ("boot-failed", with the cause), and it calls neither
    femm_session.open_femm() nor the user's initializer (user code starts processes too). Shown in a stand-in process that runs the worker's main function with femm_queue._Job
    made to raise and open_femm and an initializer that only record that they were called. The test process itself is never put into a job."""
    result = tmp_path / "result.json"
    done = subprocess.run([sys.executable, "-I", STAND_IN_POOL, "no_job", str(result)], capture_output=True, text=True, timeout=60)
    assert done.returncode == 0, done.stderr
    said = json.loads(result.read_text(encoding="utf-8"))
    assert said["kind"] == "boot-failed", said                                                # not "ready": the worker did not go on
    assert "cannot join a job of its own" in said["text"] and "cannot create a job object for the worker (a stand-in refusal)" in said["text"], said["text"]
    assert said["open_femm"] is False and said["initializer"] is False, said


def test_a_process_that_is_in_its_job_of_its_own_refuses_a_second_one_and_keeps_the_first(tmp_path):
    """The guard of _join_own_job (T13-P4: its second-call branch was never run). A process that is in its job of its own does not make another: closing the first job would end
    this very process. Shown in a stand-in process, never in the test process (a job's limits would hold pytest and every child it starts): the first call puts the stand-in into
    a job; the second raises 'already in a job of its own', and the first job is still the process's, with its handle open."""
    result = tmp_path / "result.json"
    done = subprocess.run([sys.executable, "-I", STAND_IN_POOL, "join_twice", str(result)], capture_output=True, text=True, timeout=60)
    assert done.returncode == 0, done.stderr
    said = json.loads(result.read_text(encoding="utf-8"))
    assert said["in_a_job_after_the_first"] is True, said
    assert said["second"] is not None and "already in a job of its own" in said["second"], said
    assert said["first_job_kept"] is True and said["first_job_handle_open"] is True, said


def test_a_worker_that_a_job_refuses_closes_the_job_it_made_and_joins_none_and_the_message_names_the_worker(tmp_path):
    """The other branch of _join_own_job that was never run (T13-P4): the job is made and then refuses this process (say the owner's own job allows no other). The job it made is
    closed (it holds nothing, so closing ends nobody), no job is kept, the process is in no job it was not in, and the error says 'worker' where run_solver's would say 'program'
    (T13-P5: _Job's messages reach the boot-failed text of a worker). A stand-in process again, with _Job.assign made to refuse."""
    result = tmp_path / "result.json"
    done = subprocess.run([sys.executable, "-I", STAND_IN_POOL, "assign_fails", str(result)], capture_output=True, text=True, timeout=60)
    assert done.returncode == 0, done.stderr
    said = json.loads(result.read_text(encoding="utf-8"))
    assert said["said"] is not None and "cannot join a job of its own" in said["said"] and "cannot put this worker into its job" in said["said"], said
    assert "the program" not in said["said"], said
    assert said["closed"] == ["this worker"] and said["own_job_kept"] is False, said
    assert said["in_a_job_after"] == said["in_a_job_before"], said


def test_workers_of_an_owner_that_is_itself_in_a_job_join_jobs_of_their_own_that_nest_under_it(tmp_path, lock_dir):
    """The composition question of the final review (T13-P4): the owner of a pool can be in a job already (run_solver's program is, and so is whatever a harness starts in one),
    and Windows nests jobs, so each worker's own job should become a child of the owner's and everything should go on working.

    A stand-in owner (never the test process) is put into a kill-on-close job, made here like run_solver's and never joined by this process, while it waits for a file; then it
    opens a pool of two workers (no FEMM) whose initializers start a stand-in grandchild each. Shown:
      * the pool opens: both workers joined a job of their own below the owner's (a worker that cannot says boot-failed, and the pool refuses to open);
      * every worker says it has a job of its own, is in it, and holds its handle as one that is not inheritable;
      * the outer job counts the owner, both workers and both grandchildren (the nested ones are in its accounting);
      * one worker is ended with TerminateProcess through a handle: its grandchild goes within 5 s (the worker's own job, nested, still works) and the other worker's stays;
      * the outer job is closed: the owner, the other worker and its grandchild go too (the owner's job outlives and rules the jobs below it).
    Everything the test left is ended in the finally, through a handle, by pid and creation time."""
    outer = femm_queue._Job(BELOW_NORMAL)
    launcher = subprocess.Popen([sys.executable, "-I", STAND_IN_POOL, "nested_owner", str(lock_dir), str(tmp_path)], stdin=subprocess.DEVNULL, stdout=subprocess.DEVNULL,
                                stderr=subprocess.DEVNULL)
    handles = []
    try:
        deadline = time.monotonic() + 90
        while not list(tmp_path.glob("owner_*.txt")) and time.monotonic() < deadline:
            time.sleep(0.01)
        owner = qt.Owned(*pid_and_time(next(tmp_path.glob("owner_*.txt"))))                   # the interpreter of the stand-in (the launcher in front of it is not the owner)
        handles.append(owner)
        outer.assign(owner.pid)                                                               # while it waits for 'go': nothing it starts afterwards can be born outside
        (tmp_path / "go").write_text("go", encoding="utf-8")
        deadline = time.monotonic() + 90
        while not (tmp_path / "pool_open.json").exists() and time.monotonic() < deadline:
            time.sleep(0.02)
        assert (tmp_path / "pool_open.json").exists(), "the pool did not open in 90 s: a worker could not join a job of its own below its owner's"
        opened = json.loads((tmp_path / "pool_open.json").read_text(encoding="utf-8"))
        assert len(opened["workers"]) == 2
        for pid, made, inside, inheritable in opened["facts"]:
            assert made and inside and not inheritable, f"worker {pid}: job of its own {made}, in it {inside}, handle inheritable {inheritable}"
        workers = [qt.Owned(*pid_and_time(tmp_path / f"worker_{pid}.txt")) for pid in opened["workers"]]
        handles += workers
        grandchildren = {}
        for worker in workers:
            (g_pid, g_created), = grandchildren_written_by(tmp_path, worker.pid)
            grandchildren[worker.pid] = qt.Owned(g_pid, g_created)
            handles.append(grandchildren[worker.pid])
        assert outer.active() >= 5, f"the outer job counts {outer.active()} processes: the owner, two workers and two grandchildren were expected (nested jobs are in a parent's accounting)"
        victim, survivor = workers
        victim.end()                                                                          # TerminateProcess through the handle: no clean-up at all
        assert victim.wait_gone(10), "the worker was not ended"
        assert grandchildren[victim.pid].wait_gone(5), "the grandchild of the ended worker is still running 5 s after it: a worker's job does not work below its owner's"
        assert grandchildren[survivor.pid].running() and survivor.running() and owner.running(), "something else went with the ended worker"
        outer.close()                                                                         # kill on close: the owner's job ends what is left below it
        assert owner.wait_gone(10), "the owner is still running 10 s after its job was closed"
        assert survivor.wait_gone(10) and grandchildren[survivor.pid].wait_gone(10), "a worker, or its grandchild, outlived the closing of its owner's job"
    finally:
        outer.close()                                                                         # (safe to call twice)
        launcher.kill()
        launcher.wait(30)
        for handle in handles:
            handle.end()
            handle.close()
        end_each_written_down(grandchildren_written_by(tmp_path))


def test_on_windows_the_watch_thread_only_ends_the_worker_and_leaves_the_rest_to_its_job(tmp_path):
    """The watch thread of a worker, when its owner is gone. It used to run taskkill /T /F on the worker (a walk of the tree, once); now, on Windows, it ends the worker at once
    with os._exit(1) and starts no walk: whatever the worker started is in its job, and the system ends it. Shown in a stand-in process with an owner that is gone already and an
    _end_tree that records its calls."""
    result = tmp_path / "result.json"
    done = subprocess.run([sys.executable, "-I", STAND_IN_POOL, "watch", str(result)], capture_output=True, text=True, timeout=60)
    assert done.returncode == 0, done.stderr
    said = json.loads(result.read_text(encoding="utf-8"))
    assert said["exit_code"] == 1, said
    assert said["end_tree_calls"] == [], f"the watch thread walked the worker's tree with taskkill: {said}"


def test_every_pool_worker_is_in_a_job_of_its_own_and_its_only_handle_is_not_inheritable(lock_dir):
    """Named risk: the job's handle must not reach FEMM (or any child). If a child inherited it, the job would stay open after its worker died and nothing would be ended.
    Asked of every worker, from inside: femm_queue has made its job, the worker is in it, and the handle's inherit flag is off. (The handle is made with no security
    attributes, so it is not inheritable, and subprocess closes every other handle when it starts a program.)"""
    with femm_queue.femm_pool("own jobs", workers=2, femm=False, lock_dir=lock_dir) as pool:
        pids = set(pool.pids)
        facts = pool.map(qt.own_job_facts, [0, 1, 2, 3])
    assert {pid for pid, *_ in facts} <= pids
    for pid, made, inside, inheritable in facts:
        assert made, f"worker {pid} has no job of its own"
        assert inside, f"worker {pid} is not in its own job"
        assert not inheritable, f"worker {pid} holds its job's handle as an inheritable one: a child would keep the job open"


def test_the_handle_of_a_job_is_not_inheritable():
    """The same for the job that run_solver makes (the one class both use). The job is made here and nothing is put into it: the test process is never in a job."""
    import ctypes
    from ctypes import wintypes
    k = ctypes.WinDLL("kernel32", use_last_error=True)
    k.GetHandleInformation.argtypes = [wintypes.HANDLE, ctypes.POINTER(wintypes.DWORD)]
    k.GetHandleInformation.restype = wintypes.BOOL
    job = femm_queue._Job(BELOW_NORMAL)
    try:
        flags = wintypes.DWORD()
        assert k.GetHandleInformation(job.handle, ctypes.byref(flags))
        assert not flags.value & 0x1, "the job's handle is inheritable (HANDLE_FLAG_INHERIT)"
    finally:
        job.close()


@pytest.mark.femm
def test_femm_pool_gives_every_worker_its_own_femm_below_normal_and_leaves_none_running(tmp_path, lock_dir):
    """femm=True (the default): each worker opens its own FEMM, before the user's initializer, and closes it at its end.

    Two workers, four magnet solves with private file names. The two FEMMs are two different programs, both children of their
    own worker and both at below-normal priority (the worker is lowered before FEMM starts, which FEMM inherits); the answers
    are the right ones (a magnet drives flux, no magnet drives none: nobody read the other's file); and when the pool closes, no
    FEMM is left running. This is one of the two tests here that start FEMM (the other is the next one), so it is marked `femm`.
    """
    folder = P / "out" / "tests" / "queue"
    items = [(str(folder / f"task{i}.fem"), 1_030_272 if i % 2 == 0 else 0) for i in range(4)]
    femm_pids = set()
    with femm_queue.femm_pool("real femm", workers=2, task_timeout_s=120, initializer=qt.record_femm, initargs=(str(tmp_path),),
                              lock_dir=lock_dir) as pool:
        workers = set(pool.pids)
        results = pool.map(qt.solve_on_own_femm, items)
        femm_pids = {femm_pid for _, femm_pid, _ in results}
        assert {worker for worker, _, _ in results} == workers and len(femm_pids) == 2, "each worker has a FEMM of its own"
        for worker, femm_pid, _ in results:
            assert tools.process_info(femm_pid) == ("femm_headlessw.exe", worker), "the FEMM is not the worker's own child"
            assert tools.priority_class(femm_pid) == BELOW_NORMAL
        for index, (_, _, by) in enumerate(results):
            assert (by > 0.1) if index % 2 == 0 else (abs(by) < 1e-6), (index, by)
    assert sorted(path.read_text() for path in tmp_path.glob("initializer_*.txt")) == ["femm open", "femm open"]
    assert qt.wait_until_gone(femm_pids | workers) == [], "a worker or a FEMM was left running after the pool closed"
    assert femm_queue.holder(lock_dir) is None


@pytest.mark.femm
def test_a_femm_and_its_solver_are_born_in_the_workers_job_so_the_job_can_end_them(tmp_path, lock_dir):
    """Named risk: a process that asks for CREATE_BREAKAWAY_FROM_JOB cannot start inside a job that does not allow breakaway. FEMM is started by the binding with a plain
    subprocess.Popen (no creation flags: pyfemm-qt's femm.py), and FEMM starts its own fkn.exe for each solve (closed source: not readable). So this runs the real thing: a pool
    with a real FEMM per worker, one magnet solved on it, and then the worker's own job is asked how many processes it has ever had. The solve has to have worked (a FEMM or an
    fkn.exe that could not start inside the job would fail it), and the job has to count the worker, its FEMM and the fkn.exe at least: they were born in it, which is what lets
    the job end them. Marked `femm`: it starts a FEMM."""
    folder = P / "out" / "tests" / "queue"
    with femm_queue.femm_pool("femm in a job", workers=1, task_timeout_s=120, lock_dir=lock_dir) as pool:
        ((worker, femm_pid, by, total, active),) = pool.map(qt.solve_and_count_the_job, [(str(folder / "job_task.fem"), 1_030_272)])
    assert by > 0.1, f"the magnet's flux came out as {by} T: the solve did not work inside the job"
    assert total >= 3, f"the worker's job has had {total} process(es): the worker, its FEMM and FEMM's fkn.exe were expected to be born in it"
    assert active >= 2, f"{active} process(es) in the job after the solve: the worker and its FEMM were expected"
    assert qt.wait_until_gone([worker, femm_pid]) == [], "a worker or its FEMM was left running after the pool closed"
    assert femm_queue.holder(lock_dir) is None


# ======================================================================================= run_solver
def test_run_solver_kills_on_timeout(tmp_path, lock_dir):
    """A program over timeout_s has its whole tree ended, orphans included, and run_solver returns what it had printed so far.

    The stand-in program starts a child that starts a grandchild and exits at once, so the grandchild is an orphan: the case that
    `taskkill /T` cannot reach (it finds children through their parent, and this one's has gone) and that the program's job can.
    The program and the grandchild are both gone afterwards, a bystander python.exe is not touched, and the slot has been given
    back. returncode is None: the program was ended before it had an exit code. The two processes are watched through handles
    opened as each announced itself (qt.Tracker), which also end them if the test fails.
    """
    bystander = subprocess.Popen([sys.executable, "-c", "import time; time.sleep(120)"])
    try:
        with qt.Tracker(tmp_path, ["program", "grandchild"]) as tracked:
            result = femm_queue.run_solver([sys.executable, qt.STAND_IN, "orphan", str(tmp_path)], label="stand-in solver",
                                           cwd=tmp_path, timeout_s=4, lock_dir=lock_dir)
            returncode, stdout, stderr, seconds = result                  # four values, as the plan says
            assert result.timed_out and returncode is None
            assert 4 <= seconds < 10, f"ending the job should take moments, not {seconds - 4:.1f} s after the timeout"
            assert "the program started" in stdout, "what the program wrote before it was ended should be kept"
            assert tracked.wait_announced(1), "the program and its orphaned grandchild should both have started"
            assert tracked.still_running(15) == [], "the program, or the grandchild whose parent had gone, was left running"
        assert bystander.poll() is None, "a process that run_solver did not start was ended"
        assert femm_queue.holder(lock_dir) is None
    finally:
        bystander.kill()


def test_run_solver_returns_only_when_the_whole_tree_has_ended(tmp_path, lock_dir):
    """A launcher that starts a child and exits at once (matlab.exe starting MATLAB.exe) does not free the slot early.

    The stand-in launcher exits after 0.1 s; the child it left behind holds none of its pipes (so only the job knows it is there),
    sleeps 3 s and then writes child.done. run_solver has to return after the child, with the launcher's own exit code, and say in
    `seconds` how long the tree lived. Without the job it would return in a fraction of a second and child.done would not exist.
    """
    with qt.Tracker(tmp_path, ["child"]):
        result = femm_queue.run_solver([sys.executable, qt.STAND_IN, "launcher", str(tmp_path), "3"], label="stand-in launcher",
                                       cwd=tmp_path, timeout_s=60, lock_dir=lock_dir)
        assert result.returncode == 0 and not result.timed_out, result
        assert "launcher exits now" in result.stdout
        assert (tmp_path / "child.done").read_text() == "done", "run_solver returned before the launcher's child had finished"
        assert 2.8 <= result.seconds < 30, result.seconds
        assert femm_queue.holder(lock_dir) is None


def test_a_launcher_whose_child_outlives_the_timeout_is_ended_with_it(tmp_path, lock_dir):
    """The timeout covers the whole tree: the launcher has exited, its child sleeps for ever, and after timeout_s the job ends it."""
    with qt.Tracker(tmp_path, ["child"]) as tracked:
        result = femm_queue.run_solver([sys.executable, qt.STAND_IN, "launcher", str(tmp_path), "120"], label="stand-in launcher",
                                       cwd=tmp_path, timeout_s=3, lock_dir=lock_dir)
        assert result.timed_out and result.returncode is None
        assert 3 <= result.seconds < 30, result.seconds
        assert tracked.wait_announced(1) and tracked.still_running(15) == [], "the child that outlived its launcher was left running"
        assert not (tmp_path / "child.done").exists()
        assert femm_queue.holder(lock_dir) is None


def test_a_hard_killed_caller_takes_the_program_it_started_with_it(tmp_path, lock_dir):
    """Review finding I1: end run_solver's caller with no clean-up at all, and the program and its child must not survive it.

    The caller is a process of its own that holds the slot and runs a stand-in program that has a child (both sleep). The test ends the
    caller with TerminateProcess through its handle, the way Task Manager's End task does. The operating system closes the caller's
    handles: the lock is released, and the job handle is closed, which ends the job's processes. So within a few seconds the program
    and its child are both gone, and the slot is free. Without the job the program would run on under a free slot, the next
    holder starting on top of it.
    """
    report = qt.CTX.Queue()
    caller = qt.CTX.Process(target=qt.call_run_solver, args=(str(lock_dir), str(tmp_path), "tree", 120, 900, None, report))
    try:
        with qt.Tracker(tmp_path, ["program", "child"]) as tracked:
            caller.start()
            assert tracked.wait_announced(60), "the stand-in program and its child never started"
            assert femm_queue.holder(lock_dir)["pid"] == caller.pid, "the caller should hold the slot while its program runs"
            caller.kill()                                                    # TerminateProcess: no finally, no clean-up, nothing
            caller.join(10)
            assert tracked.still_running(15) == [], "the program, or its child, outlived the caller that started it"
            assert femm_queue.holder(lock_dir) is None
        assert report.empty(), "the killed caller cannot have reported anything"
    finally:
        if caller.is_alive():
            caller.kill()


def test_ctrl_c_while_run_solvers_program_runs_leaves_no_program_and_a_free_slot(tmp_path, lock_dir):
    """A Ctrl-C while the program runs (raised the way test_ctrl_c_while_waiting_leaves_no_lock raises it) ends the whole tree.

    run_solver has to notice the interrupt (it waits in slices of half a second, because one long wait cannot be interrupted on
    Windows), end the job and wait until the job holds no process, and give the slot back, while the KeyboardInterrupt goes on to the
    caller. Afterwards the handles of the test confirm that the program and its child are gone, and the slot is free.
    (An earlier version of this test also looked at the processes the instant the interrupt reached the caller. That look raced:
    the job's count reaches zero a few microseconds before the process objects are signalled as exited, and it failed once in 25.)
    """
    report = qt.CTX.Queue()
    caller = qt.CTX.Process(target=qt.call_run_solver, args=(str(lock_dir), str(tmp_path), "tree", 120, 900, 3.0, report))
    try:
        with qt.Tracker(tmp_path, ["program", "child"]) as tracked:
            caller.start()
            assert tracked.wait_announced(60), "the stand-in program and its child never started"
            assert report.get(timeout=60) == ("interrupted",), "the Ctrl-C did not reach the caller as KeyboardInterrupt"
            caller.join(30)
            assert tracked.still_running(15) == [], "the program, or its child, was left running after the Ctrl-C"
            assert femm_queue.holder(lock_dir) is None
            with femm_queue.femm_slot("after", wait=False, lock_dir=lock_dir):
                pass
    finally:
        if caller.is_alive():
            caller.kill()


def test_a_program_that_cannot_be_put_in_its_job_is_never_run(tmp_path, lock_dir, monkeypatch):
    """If assigning to the job fails (say the caller's own job allows no other), the program is ended before it has run anything.

    run_solver raises an error that names the cause; the program (created suspended, never resumed) is dead and has written nothing;
    the slot is free. A program is never run unprotected, silently. Popen is wrapped here only to get hold of the process it makes.
    """
    made = []
    real_popen = subprocess.Popen

    class RecordingPopen(real_popen):
        def __init__(self, *args, **kwargs):
            super().__init__(*args, **kwargs)
            made.append(self)

    def refuse(self, pid):
        raise RuntimeError("cannot put the program into its job: [WinError 5] Access is denied. (simulated)")

    monkeypatch.setattr(femm_queue._Job, "assign", refuse)
    monkeypatch.setattr(subprocess, "Popen", RecordingPopen)
    code = "import pathlib, sys; pathlib.Path(sys.argv[1], 'ran.txt').write_text('it ran')"
    try:
        with pytest.raises(RuntimeError, match="into its job"):
            femm_queue.run_solver([sys.executable, "-c", code, str(tmp_path)], label="unprotected", cwd=tmp_path, timeout_s=30,
                                  lock_dir=lock_dir)
        assert len(made) == 1 and made[0].poll() is not None, "the program that could not be put in a job was left alive (suspended)"
        time.sleep(1)                                                        # a program that wrongly ran would have written by now
        assert not (tmp_path / "ran.txt").exists()
        assert femm_queue.holder(lock_dir) is None
    finally:
        for program in made:            # if run_solver failed to end it, this test must not leave a suspended program behind for ever
            if program.poll() is None:  # (it is ended through the Popen this test holds, never by a bare pid)
                program.kill()


def test_run_solver_returns_exit_code_output_and_time(tmp_path, lock_dir):
    """A program that ends by itself: its exit code, what it printed on each stream, and how long it ran. A failure is not an error."""
    code = "import sys; print('to stdout'); print('to stderr', file=sys.stderr); sys.exit(3)"
    result = femm_queue.run_solver([sys.executable, "-c", code], label="stand-in solver", cwd=tmp_path, timeout_s=60, lock_dir=lock_dir)
    assert not result.timed_out and result.returncode == 3
    assert result.stdout.strip() == "to stdout" and result.stderr.strip() == "to stderr"
    assert 0 < result.seconds < 60
    assert femm_queue.holder(lock_dir) is None


def test_run_solver_starts_the_program_below_normal_with_its_own_thread_variables(tmp_path, lock_dir):
    """The program runs at below-normal priority and sees the thread counts, and only it does: this process's environment is unchanged.

    The stand-in prints the three variables and its own priority class (what the operating system says, not what we asked).
    threads=None means the configured maximum ([femm] workers, 8 today), so a caller that forgets to choose is still capped.
    """
    code = ("import ctypes, os\n"
            "k = ctypes.WinDLL('kernel32')\n"
            "k.GetCurrentProcess.restype = ctypes.c_void_p\n"
            "k.GetPriorityClass.argtypes = [ctypes.c_void_p]\n"
            "print(os.environ.get('OMP_NUM_THREADS'), os.environ.get('CCX_NPROC_EQUATION_SOLVER'), "
            "os.environ.get('CCX_NPROC_STIFFNESS'), k.GetPriorityClass(k.GetCurrentProcess()))\n")
    names = ("OMP_NUM_THREADS", "CCX_NPROC_EQUATION_SOLVER", "CCX_NPROC_STIFFNESS")
    before = [os.environ.get(name) for name in names]
    two = femm_queue.run_solver([sys.executable, "-c", code], label="two threads", cwd=tmp_path, timeout_s=60, threads=2, lock_dir=lock_dir)
    assert two.stdout.split() == ["2", "2", "2", str(BELOW_NORMAL)], two
    default = femm_queue.run_solver([sys.executable, "-c", code], label="default threads", cwd=tmp_path, timeout_s=60, lock_dir=lock_dir)
    cap = str(paths.setting("femm", "workers"))
    assert default.stdout.split() == [cap, cap, cap, str(BELOW_NORMAL)], default
    assert [os.environ.get(name) for name in names] == before, "run_solver changed this process's own environment"


# ======================================================================================= the job holds every process at the class run_solver was asked for
# Fix round 1 of Task 13: MATLAB's graphics hosts (MATLABWindow.exe, started when a script draws a figure) were seen at normal and above normal inside run_solver's job,
# because a job that only kills its members on close does not stop a member from raising its own class, or a child from being created with a higher one. The job now carries
# JOB_OBJECT_LIMIT_PRIORITY_CLASS (0x20) beside JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE (0x2000). The model is a probe (a throw-away script, not kept in the repository, of the
# tools plan), which ran non-elevated: with the limit, SetPriorityClass(above_normal) and (normal) return 1 and the class stays, and a grandchild created with an explicit
# NORMAL_PRIORITY_CLASS flag is born at the job's class.
RAISER = r'''
import ctypes, subprocess, sys
k = ctypes.WinDLL("kernel32")
k.GetCurrentProcess.restype = ctypes.c_void_p
k.GetPriorityClass.argtypes = [ctypes.c_void_p]
k.SetPriorityClass.argtypes = [ctypes.c_void_p, ctypes.c_uint32]
NAMES = {0x20: "normal", 0x40: "idle", 0x80: "high", 0x4000: "below_normal", 0x8000: "above_normal"}
me = k.GetCurrentProcess()
print("at start:", NAMES[k.GetPriorityClass(me)], flush=True)
k.SetPriorityClass(me, 0x8000)
print("after asking for above_normal:", NAMES[k.GetPriorityClass(me)], flush=True)
k.SetPriorityClass(me, 0x20)
print("after asking for normal:", NAMES[k.GetPriorityClass(me)], flush=True)
GRANDCHILD = ("import ctypes; k = ctypes.WinDLL('kernel32'); k.GetCurrentProcess.restype = ctypes.c_void_p; k.GetPriorityClass.argtypes = [ctypes.c_void_p]; "
              "print({0x20: 'normal', 0x40: 'idle', 0x4000: 'below_normal', 0x8000: 'above_normal'}[k.GetPriorityClass(k.GetCurrentProcess())])")
done = subprocess.run([sys.executable, "-I", "-c", GRANDCHILD], creationflags=0x20, capture_output=True, text=True)     # created with an explicit NORMAL_PRIORITY_CLASS
print("grandchild created with NORMAL_PRIORITY_CLASS:", done.stdout.strip(), flush=True)
'''


def classes_the_raiser_reports(output: str) -> list[str]:
    """The class name at the end of each line the raiser printed: at start, after asking for above normal, after asking for normal, the grandchild's."""
    return [line.rsplit(":", 1)[1].strip() for line in output.strip().splitlines()]


def test_a_program_in_the_solvers_job_cannot_raise_its_own_class_or_its_childs(tmp_path, lock_dir, monkeypatch):
    """Every process of run_solver's job stays at the class it was asked for, whatever it does itself.

    The stand-in reports its class, asks the system for above normal and then for normal (both calls succeed and change nothing), and starts a grandchild with an explicit
    NORMAL_PRIORITY_CLASS flag. All four reports say below_normal. With [femm] priority = 'idle' they all say idle: the job holds the class the call was asked for, not a fixed one."""
    monkeypatch.setattr(paths, "CONFIG", write_config(tmp_path))
    below = femm_queue.run_solver([sys.executable, "-c", RAISER], label="raiser", cwd=tmp_path, timeout_s=60, lock_dir=lock_dir)
    assert below.returncode == 0, below
    assert classes_the_raiser_reports(below.stdout) == ["below_normal"] * 4, below.stdout
    monkeypatch.setattr(paths, "CONFIG", write_config(tmp_path, priority="idle"))
    idle = femm_queue.run_solver([sys.executable, "-c", RAISER], label="raiser, idle", cwd=tmp_path, timeout_s=60, lock_dir=lock_dir)
    assert idle.returncode == 0, idle
    assert classes_the_raiser_reports(idle.stdout) == ["idle"] * 4, idle.stdout
    assert femm_queue.holder(lock_dir) is None


def test_without_the_job_the_same_program_does_raise_its_class_and_its_childs_which_is_what_the_job_stops(tmp_path):
    """The control: started directly (below normal at birth, no job), the same stand-in rises to above normal when it asks, to normal when it asks again, and its grandchild is
    born at normal. So the test above can tell a job that holds the class from one that does not."""
    done = subprocess.run([sys.executable, "-c", RAISER], creationflags=subprocess.BELOW_NORMAL_PRIORITY_CLASS, capture_output=True, text=True, timeout=60)
    assert done.returncode == 0, done
    assert classes_the_raiser_reports(done.stdout) == ["below_normal", "above_normal", "normal", "normal"], done.stdout


def test_the_job_is_made_with_the_priority_class_limit_beside_the_kill_on_close_limit():
    """What the job object itself says (QueryInformationJobObject, extended limits): both limit flags, and the class it was given. `_Job` takes the class run_solver was asked for."""
    w = femm_queue._win()
    for given in (BELOW_NORMAL, IDLE):
        job = femm_queue._Job(given)
        try:
            limits = w.ExtendedLimit()
            assert w.k.QueryInformationJobObject(job.handle, 9, w.ctypes.byref(limits), w.ctypes.sizeof(limits), None)
            flags = limits.BasicLimitInformation.LimitFlags
            assert flags & 0x2000, "the job no longer ends its processes when it is closed"
            assert flags & 0x20, "the job does not hold its processes at a class"
            assert limits.BasicLimitInformation.PriorityClass == given
        finally:
            job.close()
    with pytest.raises(TypeError):
        femm_queue._Job()                                       # a job without a class would not hold one: the class is required


def test_run_solver_refuses_a_string_for_a_command_and_threads_it_may_not_use(tmp_path, lock_dir):
    """Mistakes are refused before the slot is taken or anything is started."""
    with pytest.raises(TypeError, match="list"):
        femm_queue.run_solver("ccx model", label="x", cwd=tmp_path, timeout_s=5, lock_dir=lock_dir)
    for threads in (0, 9, 2.5):
        with pytest.raises(ValueError, match="threads"):
            femm_queue.run_solver([sys.executable, "-c", "pass"], label="x", cwd=tmp_path, timeout_s=5, threads=threads, lock_dir=lock_dir)
    for timeout in (None, 0, -1):
        with pytest.raises(ValueError, match="timeout_s"):
            femm_queue.run_solver([sys.executable, "-c", "pass"], label="x", cwd=tmp_path, timeout_s=timeout, lock_dir=lock_dir)
    assert not lock_dir.exists(), "a refused call took the slot"


def test_run_solver_waits_for_a_slot_that_a_pool_or_anyone_holds(tmp_path, lock_dir, capsys):
    """While another run holds the slot, run_solver prints the waiting line and starts nothing; when it lets go, it runs."""
    marker = tmp_path / "started.txt"
    code = f"open({str(marker)!r}, 'w').write('ran'); print('ran')"
    outcome = {}

    def solver():
        outcome["result"] = femm_queue.run_solver([sys.executable, "-c", code], label="stand-in solver", cwd=tmp_path, timeout_s=60,
                                                  poll_s=0.2, lock_dir=lock_dir)

    with qt.Holder(lock_dir, "a pool, say") as holding:
        thread = threading.Thread(target=solver, daemon=True)
        thread.start()
        time.sleep(1.0)
        assert thread.is_alive() and not marker.exists(), "the program was started while another run held the slot"
        text = capsys.readouterr().out
        assert f"waiting for the solver slot: held by 'a pool, say' (pid {holding.pid})" in text
        holding.release()
        thread.join(60)
    assert not thread.is_alive() and outcome["result"].stdout.strip() == "ran" and marker.exists()
    assert femm_queue.holder(lock_dir) is None


def test_a_program_started_by_run_solver_cannot_take_the_slot_again(tmp_path, lock_dir):
    """The solver's child inherits the claim "my parent holds the slot": it is refused at once, it does not wait for its parent."""
    code = ("import pathlib, sys, femm_queue\n"
            "try:\n"
            "    with femm_queue.solver_slot('child', wait=False, lock_dir=pathlib.Path(sys.argv[1])):\n"
            "        print('TOOK THE SLOT')\n"
            "except femm_queue.QueueBusy as error:\n"
            "    print('BUSY', error)\n"
            "except RuntimeError as error:\n"
            "    print('REFUSED', error)\n")
    result = femm_queue.run_solver([sys.executable, "-c", code, str(lock_dir)], label="parent run", cwd=P, timeout_s=60, lock_dir=lock_dir)
    assert result.stdout.startswith("REFUSED"), result
    assert "holds the solver slot" in result.stdout


# ======================================================================================= the tests' own safety net
def test_a_handle_is_refused_unless_it_is_on_the_process_that_announced_itself(tmp_path):
    """Review minor 1: the tests end processes only through handles, and a handle must be on the process the test started.

    A pid names a process only while it lives; the stand-in therefore announces its creation time with it, and qt.Owned refuses a handle
    whose time (read through the handle) is not the announced one. Here a stand-in this test started (and holds with a Popen) announces
    itself; the right time is accepted, a wrong one is refused and leaves the process alone, and ending goes through the accepted handle.
    """
    sleeper = subprocess.Popen([sys.executable, qt.STAND_IN, "sleeper", str(tmp_path), "30"])
    try:
        announced = tmp_path / "program.pid"
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline and not (announced.exists() and len(announced.read_text().split()) == 2):
            time.sleep(0.05)
        pid, created = (int(part) for part in announced.read_text().split())
        with pytest.raises(ProcessLookupError, match="creation time differs"):
            qt.Owned(pid, created + 1)
        owned = qt.Owned(pid, created)
        try:
            assert owned.running() and owned.created == created
            owned.end()
            assert owned.wait_gone(10) and not owned.running()
        finally:
            owned.close()
    finally:
        sleeper.kill()


# ======================================================================================= what the command line says
def test_the_status_command_says_whether_the_slot_is_free(lock_dir, capsys):
    """`python femm_queue.py [lock folder]` prints who holds the slot (exit 0 when free, 1 when held) and takes nothing."""
    assert femm_queue.main([str(lock_dir)]) == 0
    assert "free" in capsys.readouterr().out
    with qt.Holder(lock_dir, "bench B4 ccx") as first:
        assert femm_queue.main([str(lock_dir)]) == 1
        text = capsys.readouterr().out
        assert "'bench B4 ccx'" in text and f"pid {first.pid}" in text
    assert femm_queue.main([str(lock_dir)]) == 0
