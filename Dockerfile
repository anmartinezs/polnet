# Use an official Python runtime as a parent image
FROM python:3.11-buster

SHELL [ "/bin/bash","-c" ]

# Set the working directory in the container to /app
WORKDIR /app

# Install imod 
RUN mkdir /imod_install
RUN wget https://bio3d.colorado.edu/imod/AMD64-RHEL5/imod_4.11.25_RHEL7-64_CUDA10.1.sh --no-check-certificate -O /imod_4.11.25_RHEL7-64_CUDA10.1.sh
RUN chmod +x /imod_4.11.25_RHEL7-64_CUDA10.1.sh 
RUN /imod_4.11.25_RHEL7-64_CUDA10.1.sh -y -dir /imod_install
RUN export IMOD_DIR='/imod_install/IMOD'
#RUN source /imod_install/IMOD/IMOD-linux.sh
RUN echo "#!/bin/bash" > /app/entrypoint.sh
RUN echo "source /imod_install/IMOD/IMOD-linux.sh" >> /app/entrypoint.sh
RUN echo "python ./scripts/data_gen/launch_docker.py" >> /app/entrypoint.sh
RUN chmod +xr /app/entrypoint.sh

#RUN apt-get update && apt-get install -y libx11-6 libgl1-mesa-glx
RUN apt-get update && apt-get install -y libgl1-mesa-glx libxt6

# Install any needed packages specified in requirements.txt
COPY requirements.txt /app
RUN pip install --no-cache-dir -r requirements.txt

# Copy the rest of your application's code
COPY . /app

# Set the permission for directory
RUN chmod -R 777 /app

# Run your script when the container launches
ENTRYPOINT ["/app/entrypoint.sh" ]
