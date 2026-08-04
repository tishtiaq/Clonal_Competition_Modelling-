import nbformat
from nbconvert import PythonExporter

with open("FittingModels.ipynb") as f:
    nb = nbformat.read(f, as_version=4)

exporter = PythonExporter()
source, _ = exporter.from_notebook_node(nb)

with open("FittingModels.py", "w") as f:
    f.write(source)