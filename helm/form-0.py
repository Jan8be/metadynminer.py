from traitlets import default, Unicode
from tornado import gen
from kubespawner import KubeSpawner

# c.JupyterHub.spawner_class = SillySpawner
c.MappingKernelManager.cull_idle_timeout = 259200
c.MappingKernelManager.cull_connected = False
c.MappingKernelManager.cull_busy = False
c.NotebookApp.shutdown_no_activity_timeout = 259200

c.Authenticator.allowed_users = set()
