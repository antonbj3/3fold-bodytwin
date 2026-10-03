"""Use the same CLI as the profile launcher, with bounded result recovery."""
import argparse, os, runpy, sys
p=argparse.ArgumentParser();p.add_argument('profile');p.add_argument('model');p.add_argument('--dir',required=True);p.add_argument('--title',required=True);p.add_argument('prompt')
a=p.parse_args()
assert a.model=='swarm'
os.chdir(a.dir)
sys.argv=['_run_bunny.py',a.profile,a.title,a.title]
runpy.run_path('_run_bunny.py',run_name='__main__')
