import polars as pl
import os
import threading
from bird_tool_utils import in_tempdir

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
    fifo_path = "my_fifo"
    print(os.path.abspath(fifo_path))
    os.mkfifo(fifo_path)

    # Start a new thread to write CSV data to the named pipe
    csv_writer_thread = threading.Thread(target=write_csv_to_fifo, args=(fifo_path,))
    csv_writer_thread.start()

    # Use scan_csv to parse the CSV data from the named pipe in the main thread
    df = pl.scan_csv(fifo_path)

    # Wait for the writer thread to finish
    csv_writer_thread.join()

# Print the resulting DataFrame
print(df)