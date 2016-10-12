from setuptools import setup

setup(
    name="zmat",
    author="Olav Vahtras",
    author_email="olav.vahtras@gmail.com",
    version="1.0",
    py_modules=["zmat", "molconvert"],
    scripts=["molconvert"],
    install_requires=["blocked-matrix-utils"],
    )
