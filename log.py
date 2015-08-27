#! /usr/bin/env python3

import os
import shutil
import sys
import random
import time

HOME = os.environ['HOME']
MM = os.path.join(HOME, 'mm87')
exe = 'a.out' if len(sys.argv) < 2 else sys.argv[1]

exe_path = os.path.join(MM, exe)
copied_exe_path = os.path.join(MM, 'copied_' + str(random.randint(0, 10**5)))
shutil.copy(exe_path, copied_exe_path)

def make_command(seed):
    command = "java -jar tester.jar -exec '{}' -novis -seed {}".format(copied_exe_path, seed)
    return command

def get_score(seed):
    c = make_command(seed)

    start = time.time()
    try:
        output = os.popen(c).read()
    except:
        raise "exit"

    score = float(output.split()[-1])
    if score == -1:
        score = 0
    exe_time = time.time() - start

    return {'seed': seed, 'score': score, 'time': exe_time}

def single(seeds):
    for seed in seeds:
        result = get_score(seed)
        print('{:5d} {:10.2f} {:.3f}'.format(seed, result['score'], result['time']))
        sys.stdout.flush()

def get_s(result):
    seed = result['seed']
    return '{:5d} {:10.2f} {:.3f}'.format(seed, result['score'], result['time'])

def run(seed):
    r = get_score(seed)
    sys.stderr.write(get_s(r) + '\n')
    return r
def multi(seeds):
    from multiprocessing import Pool
    pool = Pool()
    results = pool.map(run, seeds)
    for result in results:
        seed = result['seed']
        print('{:5d} {:10.2f} {:.3f}'.format(seed, result['score'], result['time']))
        sys.stdout.flush()

try:
#     single(range(1, 10))
    multi(range(1, 500))
finally:
    os.remove(copied_exe_path)

