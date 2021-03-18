from setuptools import setup, find_packages

setup(
      # mandatory
      name='banshee',
      # mandatory
      version='0.1',
      # mandatory
      author_email='joshua.bopp@stud-mail.uni-wuerzburg.de',
      packages=['banshee'],
      package_data={},
      install_requires=['pybedtools', 'click'],
      entry_points={
        'console_scripts': ['banshee = banshee.cli:start']
      }
)
