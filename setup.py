from setuptools import setup, find_packages

setup(
      # mandatory
      name='benchee',
      # mandatory
      version='0.1',
      # mandatory
      author_email='joshua.bopp@stud-mail.uni-wuerzburg.de',
      packages=['benchee'],
      package_data={},
      install_requires=['pybedtools', 'click'],
      entry_points={
        'console_scripts': ['benchee = benchee.cli:start']
      }
)
