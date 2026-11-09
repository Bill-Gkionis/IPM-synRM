"""Why FEMM workers must be separate processes, each with its own FEMM and its own files.
Run, from the learn folder: uv run python 04_workers.py

Part 1: 4 process workers -> 4 separate femm_headlessw.exe, one per worker. Each worker asks femm_session for a FEMM
        of its own, and femm_session gives every one its own free port, so nothing is shared. (With FEMM's own
        openfemm() every worker would dial port 48015, and only one of them could be connected.)
Part 2: 2 thread workers -> they fail: threads share one Python, the femm.mi_*() calls reach one FEMM per process,
        so a process can hold only one FEMM, and the second open_femm() is refused.
Part 3: 2 process workers sharing ONE file name -> worker 1 reads worker 2's answer.
        Each has its own FEMM, but both FEMMs save to and solve the same files on disk.
        The sleeps force the bad timing so it happens every run; in real runs it happens "sometimes".
Part 4: the same two jobs with private file names -> each gets its own answer.

Numbers to expect (FEMM 5.0.20261101, run on 2026-11-05, 17 s; the process ids differ on every run):
    part 1: femm_headlessw.exe running before: 0
            femm_headlessw.exe running while 4 process workers are busy: 4
            worker PIDs: [5224, 6204, 10120, 21984]
            femm_headlessw.exe running after the pool closed: 0
    part 2: thread finished, PID 14220
            thread failed: RuntimeError this process already has a FEMM open: close_femm() it first. The module-level
            femm.mi_*() calls reach one FEMM per process, so a second would leave the first one running. Threads
            share a process; give every worker its own process (see learn/04_workers.py).
            (that is one printed line, wrapped here; which of the two threads fails can change, but one always does)
    part 3: out/shared.fem: my magnet has Hc = 1,030,272 A/m, I read B = 0.0000 T
            out/shared.fem: my magnet has Hc =         0 A/m, I read B = 0.0000 T
    part 4: out/job_a.fem: my magnet has Hc = 1,030,272 A/m, I read B = 0.1662 T
            out/job_b.fem: my magnet has Hc =         0 A/m, I read B = 0.0000 T
Reading these:
  - The "4" of part 1 depends on timing: the main process counts after a 2.5 s sleep, so on a busy PC, where a worker's
    FEMM has not started by then, it can read lower. It counts femm_headlessw.exe on the whole PC, so it reads 0 / 4 / 0
    only when no other headless FEMM is running.
  - Part 3: the magnet job reads 0.0000 T although its magnet is strong, because it reads the .ans that the other job
    (no magnet) wrote into the same file name. That is the failure, made to happen every run.
  - Part 4's 0.1662 T is not learn/03's 0.1582 T for the same magnet in the same box: 03 sets an A = 0 boundary on the
    four edges of the box, so no flux crosses them, and this script sets no boundary at all.
"""
import os
import subprocess
import sys
import time
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor

# `import femm_session` must find femm_session.py, which lives one folder up. Python looks for a module in the folders
# on sys.path, and when a script is run it puts the folder OF THE SCRIPT on that list (here learn\), not the folder
# above it. So this line adds the folder above: abspath(__file__) is this file's full path, the first dirname() takes
# the file name off (learn\), the second goes up one folder (Python to FEMM), and insert(0, ...) puts it first.
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
import femm_session


def femm_count():
    """How many FEMMs are running: one femm_headlessw.exe per FEMM (it starts a short-lived fkn.exe for each solve)."""
    out = subprocess.run(["tasklist", "/FI", "IMAGENAME eq femm_headlessw.exe", "/NH"], capture_output=True, text=True).stdout
    return sum(1 for ln in out.splitlines() if ln.lower().startswith("femm_headlessw.exe"))


def hold_femm(k):
    import femm
    femm_session.open_femm()                # this worker's own FEMM, on a port of its own
    try:
        femm.newdocument(0)
        time.sleep(3)                       # stay alive so the main process can count the FEMMs
        return os.getpid()
    finally:
        femm_session.close_femm()


def magnet_job(args):
    """Draw one magnet with coercivity hc, solve, read B at its centre. delay_* force the bad timing of part 3."""
    import femm
    name, hc, delay_start, delay_read = args
    time.sleep(delay_start)
    femm_session.open_femm()                # its own FEMM, but the FILE below is whatever name the job was given
    try:
        femm.newdocument(0)
        femm.mi_probdef(0, "millimeters", "planar", 1e-8, 20, 30)
        femm.mi_getmaterial("Air")
        femm.mi_addmaterial("PM", 1.05, 1.05, hc)
        femm.mi_drawrectangle(-30, -30, 30, 30)
        femm.mi_drawrectangle(-5, -1, 5, 1)
        for (x, y), mat in (((0, 0), "PM"), ((0, 20), "Air")):
            femm.mi_addblocklabel(x, y); femm.mi_selectlabel(x, y)
            femm.mi_setblockprop(mat, 1, 0, "", 90, 0, 0); femm.mi_clearselected()
        femm.mi_saveas(os.path.abspath(name))
        femm.mi_analyze(1)                  # writes <name>.ans next to <name>.fem
        time.sleep(delay_read)              # ... meanwhile another worker may overwrite that .ans
        femm.mi_loadsolution()              # reads <name>.ans from disk
        return f"{name}: my magnet has Hc = {hc:>9,} A/m, I read B = {femm.mo_getb(0, 0)[1]:.4f} T"
    finally:
        femm_session.close_femm()


if __name__ == "__main__":                  # Windows starts each worker by re-importing this file: this guard stops
    os.makedirs("out", exist_ok=True)       # every worker from starting its own pool, forever
    print("part 1: femm_headlessw.exe running before:", femm_count())
    with ProcessPoolExecutor(4) as pool:
        futures = [pool.submit(hold_femm, k) for k in range(4)]
        time.sleep(2.5)
        print("        femm_headlessw.exe running while 4 process workers are busy:", femm_count())
        print("        worker PIDs:", sorted(f.result() for f in futures))
    time.sleep(1)
    print("        femm_headlessw.exe running after the pool closed:", femm_count())

    print("part 2: the same with 2 threads")
    with ThreadPoolExecutor(2) as pool:
        for f in [pool.submit(hold_femm, k) for k in range(2)]:
            try:
                print("        thread finished, PID", f.result())
            except Exception as e:
                print("        thread failed:", type(e).__name__, e)

    print("part 3: two process workers, one shared file name")
    jobs = [("out/shared.fem", 1_030_272, 0.0, 4.0), ("out/shared.fem", 0, 1.0, 0.0)]   # strong magnet, no magnet
    with ProcessPoolExecutor(2) as pool:
        for line in pool.map(magnet_job, jobs):
            print("       ", line)

    print("part 4: the same jobs, private file names")
    jobs = [("out/job_a.fem", 1_030_272, 0.0, 4.0), ("out/job_b.fem", 0, 1.0, 0.0)]
    with ProcessPoolExecutor(2) as pool:
        for line in pool.map(magnet_job, jobs):
            print("       ", line)
