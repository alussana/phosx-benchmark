FROM ubuntu
RUN apt-get update
RUN apt-get upgrade -y 
RUN apt-get install -y wget
RUN apt-get install -y zip
RUN apt-get install -y python3
RUN ln -s /usr/bin/python3 /usr/bin/python
RUN wget -O get-pip.py "https://bootstrap.pypa.io/get-pip.py"
RUN rm /usr/lib/python3.12/EXTERNALLY-MANAGED
RUN python get-pip.py
RUN pip install h5py
RUN pip install numpy
RUN pip install pandas
RUN pip install openpyxl
RUN pip install matplotlib
RUN pip install seaborn
RUN pip install logomaker
RUN pip install requests
RUN pip install scipy
RUN pip install scikit-learn
RUN apt-get install -y rustc
RUN apt-get install -y cargo
RUN pip install gseapy==1.0.5
RUN apt-get -y install git
COPY kinex .
RUN cd kinex & pip install .
ARG CACHEBUST=1
COPY phosx .
RUN cd phosx & pip install .

