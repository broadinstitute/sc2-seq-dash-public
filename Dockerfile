# Move to /dash

FROM gcr.io/gcid-viral-seq/plotly-dash-gcp:main

WORKDIR /

COPY . /

EXPOSE 8000

ENTRYPOINT ["gunicorn","--bind=0.0.0.0:8080","main:server"]
