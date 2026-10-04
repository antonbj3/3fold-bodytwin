"""Regression checks for idle-agent overcommit and machine-specific bounds."""
import contextlib
import importlib.util
import io
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

spec=importlib.util.spec_from_file_location('guard',Path(__file__).with_name('launch_reserved.py'))
g=importlib.util.module_from_spec(spec);spec.loader.exec_module(g)
M=g.MIB

class CapacityTests(unittest.TestCase):
    def setUp(self):
        self.temp=tempfile.TemporaryDirectory();self.root=Path(self.temp.name)
        self.slice=self.root/'slice';self.slice.mkdir()
        self.mem=self.root/'meminfo';self.mem.write_text('MemAvailable: 30000000 kB\n')
        self.limit=28*1024*M;self.current=0
        def paths(name):
            return {'/sys/fs/cgroup/research.slice':self.slice,'/proc/meminfo':self.mem,'/opt/agents/capacity_admission.lock':self.root/'lock'}[name]
        self.paths=patch.object(g,'Path',side_effect=paths);self.paths.start()
        self.systemctl=patch.object(g.subprocess,'check_output',side_effect=lambda *a,**k:f'MemoryMax={self.limit}\nMemoryCurrent={self.current}\n');self.systemctl.start()
        self.disk=patch('shutil.disk_usage',return_value=type('Disk',(),{'free':20*1024**3})());self.disk.start()
    def tearDown(self):
        self.disk.stop();self.systemctl.stop();self.paths.stop();self.temp.cleanup()
    def agent(self,name,maximum=1500*M,current=10*M):
        d=self.slice/name;d.mkdir();(d/'memory.max').write_text(str(maximum));(d/'memory.current').write_text(str(current));self.current+=current
    def test_ovh_idle_workers_reserve_full_limit(self):
        for i in range(16):self.agent(f'agent-{i}.service')
        self.assertEqual(g.snapshot(1500,2048,4096,17)['slots'],1)
        self.agent('agent-16.service')
        self.assertEqual(g.snapshot(1500,2048,4096,17)['slots'],0)
    def test_upcloud_never_oversubscribes_idle_agents(self):
        self.limit=7*1024*M
        for i in range(3):self.agent(f'agent-{i}.service')
        self.assertEqual(g.snapshot(1500,1536,1536,3)['slots'],0)
    def test_actual_larger_limit_and_non_agent_memory_count(self):
        self.agent('agent-other.service',4000*M,10*M)
        self.agent('background.service',3000*M,2000*M)
        self.current+=1000*M
        s=g.snapshot(1500,2048,4096,17)
        self.assertEqual(s['reserved'],8000*M)
        self.assertEqual(s['total'],1)
        self.assertEqual(s['slots'],12)
    def test_host_memory_and_disk_close_admission(self):
        self.mem.write_text('MemAvailable: 4000000 kB\n')
        self.assertEqual(g.snapshot(1500,2048,4096,17)['slots'],0)
        self.mem.write_text('MemAvailable: 30000000 kB\n')
        with patch('shutil.disk_usage',return_value=type('Disk',(),{'free':1})()):
            self.assertEqual(g.snapshot(1500,2048,4096,17)['slots'],0)
    def test_overbudget_blocks_launch(self):
        self.limit=7*1024*M
        for i in range(6):self.agent(f'agent-{i}.service')
        with patch('sys.argv',['guard','--host-cap','3','--','systemd-run']),patch.object(g.subprocess,'run') as run,contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(g.main(),75);run.assert_not_called()
    def test_unbounded_child_fails_closed(self):
        self.agent('agent-unbounded.service');(self.slice/'agent-unbounded.service/memory.max').write_text('max')
        with patch('sys.argv',['guard','--','systemd-run']),patch.object(g.subprocess,'run') as run,contextlib.redirect_stderr(io.StringIO()):
            self.assertEqual(g.main(),75);run.assert_not_called()
    def test_lock_is_held_until_systemd_registers_unit(self):
        import fcntl
        def launched(*a,**k):
            with (self.root/'lock').open('a') as other:
                with self.assertRaises(BlockingIOError):fcntl.flock(other,fcntl.LOCK_EX|fcntl.LOCK_NB)
            return type('Result',(),{'returncode':0})()
        with patch('sys.argv',['guard','--','systemd-run']),patch.object(g.subprocess,'run',side_effect=launched):
            self.assertEqual(g.main(),0)

if __name__=='__main__':unittest.main()
