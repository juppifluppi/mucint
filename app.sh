#!/bin/bash

streamlit run \
          --server.address 0.0.0.0 \
          --server.port 8080 \
          --server.headless True \
          mucint_streamlit.py
