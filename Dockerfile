FROM ubuntu:16.04

RUN apt-get update \
    && apt-get install -y wget
RUN wget -O /etc/apt/sources.list.d/neurodebian.sources.list http://neuro.debian.net/lists/xenial.us-ca.full
RUN apt-key adv --recv-keys --keyserver hkp://pgp.mit.edu:80 0xA5D32F012649A5A9

# Run apt-get calls
RUN apt-get update \
    && apt-get install -y fsl-5.0-core

# Install AFNI and set environment variables
RUN cd /tmp && wget --no-check-certificate https://afni.nimh.nih.gov/pub/dist/tgz/linux_openmp_64.tgz && tar xfz /tmp/linux_openmp_64.tgz && mv /tmp/linux_openmp_64 /opt/afni
ENV PATH /opt/afni:$PATH
ENV DYLD_FALLBACK_LIBRARY_PATH /opt/afni

# Configure environment
ENV FSLDIR=/usr/share/fsl/5.0
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV FSLMULTIFILEQUIT=TRUE
ENV POSSUMDIR=/usr/share/fsl/5.0
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0
ENV FSLTCLSH=/usr/bin/tclsh
ENV FSLWISH=/usr/bin/wish


# Install the MCR dependencies and some things we'll need and download the MCR
# from Mathworks - silently install it
RUN apt-get -qq update && apt-get -qq install -y unzip xorg wget curl && \
    mkdir /opt/mcr && \
    mkdir /mcr-install && \
    cd /mcr-install && \
    wget -nv http://www.mathworks.com/supportfiles/MCR_Runtime/R2012b/MCR_R2012b_glnxa64_installer.zip && \
    unzip MCR_R2012b_glnxa64_installer.zip && \
    ./install -destinationFolder /opt/mcr -agreeToLicense yes -mode silent && \
    cd / && \
    rm -rf /mcr-install

#
RUN apt-get install -y python
RUN apt-get install -y python-pip
RUN pip install nibabel

# setting relevant env variables
ENV AFNI_PATH /opt/afni
ENV FSL_PATH $FSLDIR/bin
ENV MCR_PATH /opt/mcr/v80
ENV PATH $AFNI_PATH:$FSL_PATH:$MCR_PATH:$PATH

# Configure environment variables for MCR
# ENV LD_LIBRARY_PATH /opt/mcr/v80/runtime/glnxa64:/opt/mcr/v80/bin/glnxa64:/opt/mcr/v80/sys/os/glnxa64:/opt/mcr/v80/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/opt/mcr/v80/sys/java/jre/glnxa64/jre/lib/amd64/server:/opt/mcr/v80/sys/java/jre/glnxa64/jre/lib/amd64
# ENV XAPPLRESDIR /opt/mcr/v80/X11/app-defaults

RUN mkdir -p /code
RUN mkdir /oppni
RUN mkdir /projects
RUN mkdir /scratch
RUN mkdir /local-scratch

COPY compiled/run_oppni.sh /oppni/
COPY compiled/oppni /oppni/
COPY bids/oppni.py /oppni/
ENV OPPNI /oppni
ENV PATH $OPPNI:$PATH

ENTRYPOINT ["python", "/oppni/oppni.py"]

