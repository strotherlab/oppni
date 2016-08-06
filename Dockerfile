FROM ubuntu:trusty

# Configure environment
ENV FSLDIR=/usr/share/fsl/5.0
ENV FSLOUTPUTTYPE=NIFTI_GZ
ENV FSLMULTIFILEQUIT=TRUE
ENV POSSUMDIR=/usr/share/fsl/5.0
ENV LD_LIBRARY_PATH=/usr/lib/fsl/5.0
ENV FSLTCLSH=/usr/bin/tclsh
ENV FSLWISH=/usr/bin/wish

ENV AFNI_PATH /opt/afni/
ENV FSL_PATH $FSLDIR/bin/
ENV MCR_PATH /opt/mcr/v80/
ENV PATH $AFNI_PATH:$FSL_PATH:$MCR_PATH:$PATH

ENV PATH /opt/afni:$PATH
ENV DYLD_FALLBACK_LIBRARY_PATH /opt/afni

COPY tmp/cpac_install.sh /tmp/cpac_install.sh
RUN /tmp/cpac_install.sh -s
RUN /tmp/cpac_install.sh -p 
RUN /tmp/cpac_install.sh -n afni
RUN /tmp/cpac_install.sh -n fsl


#
RUN apt-get install -y python


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

