from setuptools import setup, find_packages

setup(name='ddot',
      version='0.1',
      description='Data-Driven Ontology Toolkit',
      url='http://github.com/michaelkyu/ontology',
      author='Michael Ku Yu',
      author_email='michaelyu@alum.mit.edu',
      license='MIT',
      classifiers=[
          # How mature is this project? Common values are 
          #   3 - Alpha
          #   4 - Beta
          #   5 - Production/Stable
          'Development Status :: 3 - Alpha',

          # Indicate who your project is intended for
          'Intended Audience :: Developers',
          'Topic :: Software Development :: Build Tools',

          # Pick your license as you wish (should match "license" above)
          'License :: OSI Approved :: MIT License',

          # Specify the Python versions you support here. In particular, ensure
          # that you indicate whether you support Python 2, Python 3 or both.
          'Programming Language :: Python :: 2.7'
      ],
      keywords='ontology hierarchy',
      packages=['ddot'],
      install_requires=['pandas>=0.20', 'numpy', 'scipy', 'ndex-dev==3.0.11.35', 'python-igraph', 'networkx>=2.0', 'requests', 'tulip-python'],
      include_package_data=True,
      zip_safe=False)

      # package_data={
      #     '':['config.py', 'site.cfg'],
      # },
