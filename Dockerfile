ARG FROM_IMAGE=gitlab-registry.cern.ch/batch-team/dask-lxplus/lxdask-al9:latest
FROM ${FROM_IMAGE}

ARG CLUSTER=lxplus

ADD . . 

RUN echo "=======================================" && \
    echo "Installing HiggsDNA" && \
    echo "=======================================" && \
    if [[ ${CLUSTER} == "lxplus" ]]; then \
    yum -y update && \
    yum -y install git-lfs && \
    echo "Fixing dependencies in the image" && \
    conda install -y numba>=0.57.0 llvmlite==0.40.0 numpy>=1.22.0 && \
    pip install --upgrade dask-lxplus && \
    installed_version=$(pip show pyarrow | grep 'Version:' | awk '{print $2}'); \
    if [ "$(printf '%s\n' "$installed_version" "11.0.0" | sort -V | head -n1)" = "11.0.0" ]; then \
        pip install --upgrade pyarrow; \
    fi; \
    fi && \
    pip3 install .
