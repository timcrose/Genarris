



from ibslib.calculators import Slurm,tin_arguments

# Specify file_path for Slurm file to be submitted to queue
file_path = "test.sh"
# Specify if the program can overwrite submission scripts
overwrite = False
# Use default arguments for slurm
slurm_arguments = tin_arguments
# Must include a command for the slurm script
slurm_arguments["command"] = "echo HELLO_WORLD"

# Initialize
slurm = Slurm(tin_arguments)
# Write file
slurm.write(file_path)
# Submit the calculation
slurm.submit(file_path)

# Or write and submit the file together using calc
slurm.calc(file_path)
