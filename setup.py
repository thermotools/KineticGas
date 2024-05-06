# from setuptools import setup
# from pathlib import Path

# this_dir = Path(__file__).parent
# readme = (this_dir / 'pykingas/README.md').read_text()
# setup(
#     name='pykingas',
#     version='2.0.0',
#     packages=['pykingas'],
#     package_data={'pykingas': ['KineticGas*', 'fluids/*']},
#     description='Revised Enskog theory for Mie fluids, and other spherical potentials. Allows prediction of transport '
#                 'coefficients such as diffusion coefficients, viscosities, thermal diffusion coefficients'
#                 ' and thermal conductivities. In dense, multicomponent gas mixtures and supercritical mixtures.',
#     long_description=readme,
#     long_description_content_type='text/markdown',
#     author='Vegard Gjeldvik Jervell',
#     author_email='vegard.g.j@icloud.com',
#     url='https://github.com/thermotools/KineticGas',
#     install_requires=["numpy~=1.22",
#                       "scipy~=1.7",
#                       "thermopack~=2.1"]
# )