FROM node:22-bookworm-slim AS frontend

WORKDIR /app/web_app/frontend
COPY web_app/frontend/package*.json ./
RUN npm ci
COPY web_app/frontend ./
RUN npm run build


FROM python:3.11-slim

ENV PYTHONDONTWRITEBYTECODE=1
ENV PYTHONUNBUFFERED=1
ENV RMC_TOOLKITS_DATA_ROOT=/app

WORKDIR /app

COPY web_app/backend/requirements.txt /tmp/requirements.txt
RUN pip install --no-cache-dir -r /tmp/requirements.txt

COPY rmc_toolkits ./rmc_toolkits
COPY web_app/backend ./web_app/backend
COPY data ./data
COPY assets ./assets
COPY README.md ./README.md
COPY --from=frontend /app/web_app/frontend/dist ./web_app/frontend/dist

EXPOSE 5000

CMD ["sh", "-c", "gunicorn --bind 0.0.0.0:${PORT:-5000} --timeout 120 web_app.backend.app:app"]
