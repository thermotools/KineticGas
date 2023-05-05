import json
import os

ROOT = os.path.dirname(__file__) + '/../'

files = sorted(os.listdir(ROOT + 'pykingas/fluids'))

ostr = '# Fluid identifiers\n\n'
ostr += '| Fluid name | Fluid identifyer | CAS |\n'
ostr += '| ---------- | ---------------- | --- |\n'
for file in files:
    comp_id = file.split('.')[0]
    data = json.load(open(ROOT + 'pykingas/fluids/' + file, 'r'))
    name = data['name'][0] + data['name'][1:].lower()
    ostr += f"| {name} | {comp_id} | {data['cas_number']} |\n"

with open(ROOT + 'docs/markdown/fluid_identifiers.md', 'w') as file:
    file.write(ostr)