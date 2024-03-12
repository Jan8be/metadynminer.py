from traitlets import default, Unicode
from tornado import gen
from kubespawner import KubeSpawner

class SillySpawner(KubeSpawner):
  form_template = Unicode("""
  <h3>Home</h3>
  <div id="phomeDiv" style="display: block;">
  <input type="checkbox" id="delhome" name="delhome" value="delete">
  <label for="delhome">Erase if home exists</label><br/>
  <p style="background-color:orange;">Take care of checking this button, it removes whole home directory and previous data will be lost. Use in case only when notebook is broken so it does not start, in other cases, remove data from terminal.</p>
  </div>
  """,)

  option_template = Unicode("""
      <option value="{item}">{item}</option>""",
      config = True, help = "Template for html form options."
  )

    
  async def get_options_form(self):
    return self.form_template


c.JupyterHub.spawner_class = SillySpawner
c.MappingKernelManager.cull_idle_timeout = 259200
c.MappingKernelManager.cull_connected = False
c.MappingKernelManager.cull_busy = False
c.NotebookApp.shutdown_no_activity_timeout = 259200
