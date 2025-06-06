# built-in
import typing


def fasta_file_parser(fasta: str) -> dict[str, str]:
    with open(fasta, "r") as file:
        temp = [x.strip() for x in file.readlines()]
    fast_index = [i for i, j in enumerate(temp) if ">" in j]
    if len(fast_index) == 1:
        return {temp[0]: "".join(temp[1:])}
    if len(fast_index) >= 2:
        fasta_dict = {}
        for i, temp_index in enumerate(fast_index[:-1]):
            fasta_dict[temp[temp_index].strip(">")] = "".join(
                temp[temp_index + 1 : fast_index[i + 1]]
            )
        fasta_dict[temp[fast_index[-1]].strip(">")] = "".join(
            temp[fast_index[-1] + 1 :]
        )
        return fasta_dict
    raise RuntimeError("Could not parse file due to lack of fasta headers")


def fasta_parser(fasta: typing.TextIO) -> dict[str, str]:
    temp = [x.strip() for x in fasta]
    fast_index = [i for i, j in enumerate(temp) if ">" in j]
    if len(fast_index) == 1:
        return {temp[0]: "".join(temp[1:])}
    if len(fast_index) >= 2:
        fasta_dict = {}
        for i, temp_index in enumerate(fast_index[:-1]):
            fasta_dict[temp[temp_index].strip(">")] = "".join(
                temp[temp_index + 1 : fast_index[i + 1]]
            )
        fasta_dict[temp[fast_index[-1]].strip(">")] = "".join(
            temp[fast_index[-1] + 1 :]
        )
        return fasta_dict
    raise RuntimeError("Could not parse file due to lack of fasta headers")
