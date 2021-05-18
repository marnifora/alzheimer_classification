import tracemalloc
import subprocess
import time
import subprocess
import resource

command = "python boruta_classification.py -class -borutarun 1 -dataset rosmap /mnt/chr11/Data/rosmap/ " \
          "-num_cores 1 -perc 90"
print('COMMAND: {}'.format(command))
tracemalloc.start()
t = time.time()
try:
    output = subprocess.check_output(
                command,
                shell=True,
                stderr=subprocess.STDOUT,
                universal_newlines=True,
            )
except subprocess.CalledProcessError as exc:
    print('analysis script failed with error code {} and output:\n{}'.format(exc.returncode, exc.output))
else:
    print('analysis script succeeded with output:\n{}'.format(output))
print('TIME: {} min'.format((time.time()-t)/60))
print(resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss)
current, peak = tracemalloc.get_traced_memory()
print(f"MEMORY usage is {current / 10**6}MB; Peak was {peak / 10**6}MB")
tracemalloc.stop()
