#!/bin/bash

#set -x

source ./run_common.sh

example_name=introduction_ex3

message_running "$example_name" 
run_example "$example_name"
message_done_running "$example_name"
