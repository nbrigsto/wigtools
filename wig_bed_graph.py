import numpy as np
import matplotlib.pyplot as plt
import sys

def normalize_by_max(array):
    if array.size == 0:
        return array
    max_val = np.max(array)
    if max_val == 0:
        return array
    return array / max_val

def aggregate_tss_coverage(bed_file, wig_file, window_size):
    master = np.zeros(window_size * 2 + 1)
    temp = np.zeros(window_size * 2 + 1)
    counts = np.zeros(window_size * 2 + 1)

    with open(bed_file, 'r') as bf, open(wig_file, 'r') as wf:
        wig_line = wf.readline().strip()
        bed_line = bf.readline().strip()
        

        wig_chrom = None
        wig_pos = None
        wig_step = None
        wig_span = None

        while bed_line:
            fields = bed_line.split()
            bed_chrom = fields[0]
            bed_start = int(fields[1])
            direction = fields[5]

            while wig_line:
                if wig_line.startswith('fixedStep') or wig_line.startswith('variableStep'):
                    header = wig_line.split()

                    prev_chrom = wig_chrom 
                    wig_pos = None
                    wig_step = None
                    wig_span = None

                    for entry in header:
                        if entry != "fixedStep" and entry != "variableStep":
                            precursor = entry.split('=')[0]
                            if precursor == "chrom":
                                wig_chrom = entry.split('=')[1]
                            elif precursor == "start":
                                wig_pos = int(entry.split('=')[1])
                            elif precursor == "step":
                                wig_step = int(entry.split('=')[1])
                            elif precursor == "span":
                                wig_span = int(entry.split('=')[1])

                    print(f"We just ran into chromosome {wig_chrom}")

                    if bed_chrom == prev_chrom:
                        while bed_chrom == prev_chrom:
                            bed_line = bf.readline().strip()
                            fields = bed_line.split()
                            if len(fields) >= 6:
                                bed_chrom = fields[0]
                                bed_start = int(fields[1])
                                direction = fields[5]
                            else:
                                return master

                                    
                elif wig_chrom == bed_chrom:
                    if wig_step:
                        if wig_pos > bed_start - window_size - wig_step and wig_pos <= bed_start + window_size:
                            jmp_back_pt = wf.tell()
                            jmp_back_pos = wig_pos

                            while wig_pos <= bed_start + window_size and wig_line and not wig_line.startswith("fixedStep"):
                                index = wig_pos - bed_start + window_size
                                if wig_span:
                                    for i in range(0, wig_span):
                                        if 0 <= index + i <= window_size * 2:
                                            if direction == "+":
                                                temp[index + i] = float(wig_line)
                                                counts[index + i] += 1
                                            else:
                                                temp[window_size * 2 - index - i] = float(wig_line)
                                                counts[window_size * 2 - index + i] += 1
                                else:
                                    for i in range(0, wig_step):
                                        if 0 <= index + i <= window_size * 2:
                                            if direction == "+":
                                                temp[index + i] = float(wig_line)
                                                counts[index + i] += 1
                                            else:
                                                temp[window_size * 2 - index - i] = float(wig_line)
                                                counts[window_size * 2 - index + i] += 1

                                wig_line = wf.readline().strip()
                                wig_pos += wig_step

                            wf.seek(jmp_back_pt)
                            wig_pos = jmp_back_pos

                            #temp = normalize_by_max(temp)
                            master += temp
                            temp = np.zeros(window_size * 2 + 1)
                            break

                        elif wig_pos > bed_start + window_size:
                            wig_pos += wig_step
                            break

                        wig_pos += wig_step

                    else:

                        fields = wig_line.split()

                        wig_pos = int(fields[0])
                        wig_val = float(fields[1])

                        if not wig_span:
                            wig_span = 1

                        if wig_pos > bed_start - window_size - wig_span and wig_pos <= bed_start + window_size:
                            jmp_back_pt = wf.tell()
                            jmp_back_pos = wig_pos
                                                                            
                            while wig_pos <= bed_start + window_size and not wig_line.startswith("variableStep"):
                                index = wig_pos - bed_start + window_size
                                
                                for i in range(0, wig_span):
                                    if 0 <= index + i <= window_size * 2:
                                        if direction == "+":
                                            temp[index + i] = wig_val 
                                            counts[index + i] += 1 
                                        else:
                                            temp[window_size * 2 - index - i] = wig_val  # Sum values instead of overwriting
                                            counts[window_size * 2 - index - i] += 1

                                wig_line = wf.readline().strip()
                                if not wig_line:
                                    break
                                
                                fields = wig_line.split()
                                if len(fields) < 2:
                                    break  # Handle malformed lines
                                wig_pos = int(fields[0])
                                wig_val = float(fields[1])
                            
                            wf.seek(jmp_back_pt)
                            wig_pos = jmp_back_pos

                            #temp = normalize_by_max(temp)
                            master += temp
                            temp = np.zeros(window_size * 2 + 1)
                            break

                wig_line = wf.readline().strip()

            bed_line = bf.readline().strip()
    
    return master, counts

def plot_aggregate_coverage(master_list, window, output_file="aggregate_coverage.png"):
    """Plot the aggregate TSS coverage and save it as a .png."""
    x = np.arange(-window, window + 1)
    plt.plot(x, master_list)
    plt.xlabel("Position relative to TSS (bp)")
    plt.ylabel("Aggregate normalized coverage")
    plt.title("Aggregate TSS Coverage")
    plt.grid()
    plt.savefig(output_file, format='png', dpi=300)  # Save the plot as a .png
    plt.show()  # Display the plot


if len(sys.argv) != 6:
    print("Usage: python script.py <tss_file_path> <wig_file_path> <window_size> <png_name.png> <type>")
    sys.exit(1)

# Take command-line inputs
tss_file = sys.argv[1]
wig_file = sys.argv[2]

try:
    window_size = int(sys.argv[3])
except ValueError:
    print("Error: window_size must be an integer.")
    sys.exit(1)

# Compute aggregate coverage while iterating through the .wig file
master_list, counts_list = aggregate_tss_coverage(tss_file, wig_file, window_size)

for i in range(0, window_size * 2 + 1):
    master_list[i] = master_list[i] / counts_list[i]

# Plot the result
plot_aggregate_coverage(master_list, window_size, sys.argv[4])
with open(f"arrays_{sys.argv[5]}.txt", "a") as file:
    for item in master_list:
        file.write(f"{item}\t")
    file.write("\n")
