#!/usr/bin/env bash

wget -c -r ftp://ftp.ccb.jhu.edu/pub/RNAseq_protocol
ln -s ftp.ccb.jhu.edu/pub/RNAseq_protocol/chrX_data.tar.gz
tar xvzf chrX_data.tar.gz

