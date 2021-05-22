from getpass import getpass
import os

import click
from flask.cli import FlaskGroup
from passlib.apache import HtpasswdFile
import pkg_resources

from .flask import app

@app.cli.command()
@click.argument('username', required=False)
def passwd(username):
    """Set the password for a user."""
    if username is None:
        username = input('Username: ')
    password = getpass()

    path = os.path.join(app.instance_path, 'htpasswd')
    os.makedirs(app.instance_path, exist_ok=True)
    try:
        htpasswd = HtpasswdFile(path)
    except FileNotFoundError:
        htpasswd = HtpasswdFile()

    htpasswd.set_password(username, password)
    htpasswd.save(path)
