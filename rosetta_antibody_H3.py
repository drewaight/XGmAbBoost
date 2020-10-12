import progressbar
import fileinput  
import timeit
import subprocess
    
bar = progressbar.ProgressBar(max_value=progressbar.UnknownLength).start()
start_time = timeit.default_timer()
p = subprocess.Popen(['mpiexec', '-n', '16', 'antibody_H3.mpi.linuxgccrelease', '@flags'], 
                stdout=subprocess.PIPE, stderr=subprocess.PIPE)
result = []
while p.stdout is not None:
    bar.update()
    line = p.stdout.readline()
    result.append(line.decode('UTF-8').rstrip('\r'))
    if not line:
        print("\n")
        p.stdout.flush()
        break  
with open("antibody_H3.log", 'w') as f:
    f.write(''.join(result))
stop_time = timeit.default_timer()
total_time = stop_time - start_time
with open("time.txt", "w") as f:
    f.write(total_time)