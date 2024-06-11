FROM ubuntu:latest

RUN apt-get update && \
    apt-get install -y \
    python3 \
    python3-pip \
    python3-venv \
    openjdk-11-jre-headless \
    wget && \
    apt-get clean

COPY requirements.txt /home/ubuntu/requirements.txt

RUN wget -qO- https://get.nextflow.io | bash && \
    mv nextflow /usr/local/bin/ && \
    chmod +x /usr/local/bin/nextflow

RUN python3 -m venv /opt/venv && \
    /opt/venv/bin/pip install --upgrade pip && \
    /opt/venv/bin/pip install -r /home/ubuntu/requirements.txt

#the following part downloads and unzips the cuttleflow code, could be removed
RUN wget https://github.com/Nine-s/CuttleFlow/archive/refs/heads/master.zip -O /home/ubuntu/Cuttleflow.zip
RUN apt install unzip
RUN unzip /home/ubuntu/Cuttleflow.zip -d /home/ubuntu/Cuttleflow

ENV PATH="/opt/venv/bin:$PATH"

RUN echo 'alias python="python3"' >> /root/.bashrc
RUN echo 'alias pip="pip3"' >> root/.bashrc
CMD ["bash"]

