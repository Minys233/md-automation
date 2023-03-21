import unittest
import logging
from subprocess_runner import SubprocessRunner

logging.basicConfig(format='%(asctime)s-%(levelname)s-%(message)s', datefmt='%d-%b-%y %H:%M:%S', level=logging.INFO)

class TestSubprocessRunner(unittest.TestCase):
    def test_append_cmd(self):
        sr = SubprocessRunner('echo -n Hello')
        sr.append_cmd(['World'])
        sr.append_cmd(' YO!')
        self.assertEqual(sr.cmd, ['echo', '-n', 'Hello', 'World', 'YO!'])

    def test_run(self):
        sr = SubprocessRunner('echo -n Hello')
        sr.append_cmd(['World', 'YO!'])
        sr.run()
        self.assertEqual(sr.stdout, 'Hello World YO!')

    def test_operated(self):
        sr = SubprocessRunner('echo -n Hello')
        sr.append_cmd(['World', 'YO!'])
        self.assertEqual(sr.operated, False)
        sr.run()
        self.assertEqual(sr.operated, True)
        sr.run(ignore_opreated=True)
        self.assertEqual(sr.operated, True)

    def test_cmd(self):
        sr = SubprocessRunner('echo -n Hello')
        self.assertEqual(sr.cmd, ['echo', '-n', 'Hello'])

    def test_other_kwargs(self):
        sr = SubprocessRunner('echo -n Hello', check=True, shell=True, capture_output=True)
        sr.run()
        self.assertEqual(sr.success, True)
        self.assertEqual(sr.check, True)
        self.assertEqual(sr.stdout, "Hello")
        print(sr.proc)

if __name__ == '__main__':
    unittest.main()