FROM debian
RUN apt-get update && \
        apt-get upgrade -y && \
        apt-get install -y wget zip python3-dev
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN wget -O get-pip.py "https://bootstrap.pypa.io/get-pip.py"
RUN python get-pip.py --break-system-packages
RUN pip install --break-system-packages h5py
RUN pip install --break-system-packages "numpy<2"
RUN pip install --break-system-packages pandas
RUN pip install --break-system-packages openpyxl
RUN pip install --break-system-packages matplotlib
RUN pip install --break-system-packages seaborn
RUN pip install --break-system-packages logomaker
RUN pip install --break-system-packages requests
RUN pip install --break-system-packages scipy
RUN pip install --break-system-packages scikit-learn
RUN apt-get install -y rustc
RUN apt-get install -y cargo
RUN pip install --break-system-packages gseapy==1.0.5
RUN apt-get -y install git
COPY kinex .
RUN cd kinex & pip install --break-system-packages .
COPY phosx .
RUN cd phosx & pip install --break-system-packages .
RUN apt-get install -y r-base
RUN pip install --break-system-packages kstar
RUN wget -O KSTAR_graphs.tar.gz https://figshare.com/ndownloader/files/28768155 && \
        tar -xzf KSTAR_graphs.tar.gz # this creates /NETWORKS/NetworKIN necessary for kstar.config.update_network_directory('/NETWORKS/NetworKIN') in run_kstar.py
COPY kstar_config.py .
RUN python kstar_config.py
RUN chmod -R 755 /usr/local/lib/python3.11/dist-packages/RESOURCE_FILES.tar.gz
RUN chmod -R 755 /usr/local/lib/python3.11/dist-packages/RESOURCE_FILES
RUN chmod -R 755 /usr/local/lib/python3.11/dist-packages/kstar
