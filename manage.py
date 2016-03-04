#!/usr/bin/env python2

from flask.ext.script import Manager
from uclex import app
import uclex
#import mongodb

manager = Manager(app)


@manager.command
def load_db(): uclex.load_db()


@manager.command
def load_base_coverage(): uclex.load_base_coverage()


@manager.command
def load_variants_file(): uclex.load_variants_file()


@manager.command
def load_gene_models(): uclex.load_gene_models()


@manager.command
def load_dbsnp_file(): uclex.load_dbsnp_file()


@manager.command
def load_constraint_information(): uclex.load_constraint_information()


@manager.command
def load_mnps(): uclex.load_mnps()


@manager.command
def create_cache(): uclex.create_cache()


@manager.command
def precalculate_metrics(): uclex.precalculate_metrics()

@manager.command
def mrc_hpo(): uclex.mrc_hpo()

if __name__ == "__main__":
    manager.run()

