from setuptools import setup, find_packages

#, read_configuration

# import ConfigParser
# Config = ConfigParser.ConfigParser()
# Config.read('ddot/site.cfg')
# print 'SECTIONS', Config.sections()
# print Config.items('executables')

#conf_dict = read_configuration('ddot/setup.cfg')
#print conf_dict

#import ddot.config

# with open('tmp.py', 'w') as f:
#     f.write('x = 5\n')

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
      install_requires=['pandas', 'numpy', 'scipy', 'ndex-dev', 'python-igraph', 'networkx==1.11', 'cxmate', 'requests'],
      include_package_data=True,
      zip_safe=False)

      # package_data={
      #     '':['config.py', 'site.cfg'],
      # },
