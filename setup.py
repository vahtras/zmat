from setuptools import setup

setup(
    name="zmat",
    author="Olav Vahtras",
    author_email="olav.vahtras@gmail.com",
    version="1.0",
    install_requires=["util"],
    dependency_links=["git+https://github.com/vahtras/util.git@master#egg=util"],
    )
