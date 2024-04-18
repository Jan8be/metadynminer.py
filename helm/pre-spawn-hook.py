import asyncio 
import kubernetes_asyncio
from kubernetes_asyncio import config, client
import datetime

ns = 'metadynminer-ns'
usagelog = '/srv/jupyterhub/usage.log'

async def bootstrap_pre_spawn(spawner):
  config.load_incluster_config()
  namespace = spawner.namespace
  username = spawner.user.name
  original = username
  if "-" in username:
      username = username.replace("-", "-2d")
  if "_" in username:
      username = username.replace("_", "-5f")

  with open(usagelog,"a") as l:
    l.write(datetime.datetime.now(datetime.timezone.utc).strftime('%c %z:') + username + '\n')

  spawner.image = get_config('hub.config.notebookImage')
  spawner.cpu_limit = 1.
  spawner.cpu_guarantee = .2
  spawner.mem_limit = '8G'
  spawner.mem_guarantee = '4G'
  spawner.container_security_context = {"capabilities": {"drop": ["ALL"]}}

c.KubeSpawner.pre_spawn_hook = bootstrap_pre_spawn
c.KubeSpawner.automount_service_account_token = False
