#!/usr/bin/env python3

from kstar import config

config.install_resource_files()
config.NETWORK_DIR, config.NETWORK_Y_PICKLES, config.NETWORK_ST_PICKLES = config.update_network_directory('/NETWORKS/NetworKIN')
