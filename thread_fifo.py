import polars as pl
import os
import sys
import threading
from multiprocessing import Process
from bird_tool_utils import in_tempdir

type = sys.argv[1]
split = sys.argv[2]
"""
Testing:
python thread_fifo.py read_csv process => BrokenPipeError/no such device
python thread_fifo.py read_csv thread => BrokenPipeError/no such device
python thread_fifo.py scan_csv process => OSError: Illegal seek (os error 29)
python thread_fifo.py scan_csv thread => hangs
python thread_fifo.py base process => works
python thread_fifo.py base thread => works
"""

# Define a generator function that yields CSV strings
def csv_data_generator():
    yield "id,name\n"
    yield "1,Alice\n"
    yield "2,Bob\n"
    yield "3,Charlie\n"
    yield "4,David\n"

# Define a function to write CSV data from the generator to the FIFO
def write_csv_to_fifo(fifo_path):
    with open(fifo_path, "w") as fifo:
        for csv_string in csv_data_generator():
            fifo.write(csv_string)

with in_tempdir():
    # Create a named pipe using mkfifo
    fifo_path = os.path.abspath("my_fifo")
    print(fifo_path)
    os.mkfifo(fifo_path)

    # Start a new thread to write CSV data to the named pipe
    if split == "process":
        csv_writer_thread = Process(target=write_csv_to_fifo, args=(fifo_path,))
    else:
        csv_writer_thread = threading.Thread(target=write_csv_to_fifo, args=(fifo_path,))
    csv_writer_thread.start()

    if type == "scan_csv":
        # Use scan_csv to parse the CSV data from the named pipe in the main thread
        df = pl.scan_csv(fifo_path).collect()
    elif type == "read_csv":
        df = pl.read_csv(fifo_path)
    else:
        df = []
        with open(fifo_path, "r") as f:
            for line in f:
                df.append(line)

    # Wait for the writer thread to finish
    csv_writer_thread.join()

# Print the resulting DataFrame
print(df)
