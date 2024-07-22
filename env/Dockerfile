FROM debian
RUN apt-get update
RUN apt-get upgrade -y 
RUN apt-get install -y wget
RUN apt-get install -y zip
RUN apt-get install -y python3
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
ARG CACHEBUST=1
COPY phosx .
RUN cd phosx & pip install --break-system-packages .
RUN apt-get install -y r-base
