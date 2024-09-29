#!/bin/bash

# Check if the number of iterations is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 n"
    exit 1
fi

n=$1

# Create the messages directory if it doesn't exist
mkdir -p ./messages

# Run the Python script n times and overwrite the output files
for ((i=0; i<n; i++)); do
    python3 generate_message.py > "./messages/message$i.txt"
done

echo "Messages generated in ./messages/ directory."

