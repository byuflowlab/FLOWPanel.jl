"""Exercise the installed perf FIFO protocol before expensive fixture assembly."""
import os
import select
import sys

control = os.open(sys.argv[1], os.O_RDWR)
ack = os.open(sys.argv[2], os.O_RDWR)


def command(value):
    os.write(control, (value + "\n").encode())
    if not select.select([ack], [], [], 15)[0]:
        raise RuntimeError(f"perf acknowledgement timed out: {value}")
    if os.read(ack, 4) != b"ack\n":
        raise RuntimeError("invalid perf acknowledgement")


try:
    command("disable")
    command("enable")
    result = sum(i * i for i in range(1_000_000))
    command("disable")
    assert result > 0
finally:
    os.close(control)
    os.close(ack)
print("PASS perf enable/disable acknowledgements")
