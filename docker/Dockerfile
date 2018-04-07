# Use an official Python runtime as a parent image
FROM continuumio/miniconda3

# Set the working directory to /app
WORKDIR /app

# Copy the current directory contents into the container at /app
ADD . /app

# Install any needed packages specified in requirements.txt
RUN conda update conda
RUN conda install -c conda-forge python-igraph
RUN conda install pandas numpy scipy networkx=1.11
RUN pip install tulip-python ndex-dev

# Make port 80 available to the world outside this container
#EXPOSE 80

# Define environment variable
ENV NAME World

# Run app.py when the container launches
#CMD ["python", "app.py"]