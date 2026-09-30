FROM mambaorg/micromamba:2-ubuntu24.04
LABEL description="nf-diann container 0.3, DIA-NN 2.6.1 and pyArrow 23.0"

ARG NEW_MAMBA_USER=mambauser
ARG NEW_MAMBA_USER_ID=1
ARG NEW_MAMBA_USER_GID=1
USER root

RUN apt update && apt upgrade -y
# to have envsubst, ps

RUN apt install -y gettext-base procps
# to get DIA-NN
RUN apt install -y wget unzip libgomp1 locales dotnet-sdk-8.0

# for dinosaur (need openjdk < 25)
RUN micromamba install -y -n base -c conda-forge -c bioconda pyarrow=23.0

# For DIA-NN
# Configure locale to avoid runtime errors
RUN sed -i 's/# en_US.UTF-8/en_US.UTF-8/' /etc/locale.gen && \
    locale-gen && \
    update-locale LANG=en_US.UTF-8

# Set environment variables for locale
ENV LANG=en_US.UTF-8
ENV LANGUAGE=en_US:en
ENV LC_ALL=en_US.UTF-8

# Download DIA-NN version <version>
RUN wget https://github.com/vdemichev/DiaNN/releases/download/2.0/DIA-NN-2.6.1-Academia-Linux.zip -O /tmp/diann.Linux.zip
# Unzip the DIA-NN package
RUN cd / && unzip /tmp/diann.Linux.zip
RUN rm /tmp/diann.Linux.zip

# Set appropriate permissions for the DIA-NN folder
RUN chmod -R 775 /diann-2.6.1
ENV PATH="${PATH}:/diann-2.6.1:/opt/conda/bin"
