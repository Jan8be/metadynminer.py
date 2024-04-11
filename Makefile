#image?=cerit.io/ljocha/metadynminer
image?=ghcr.io/jan8be/metadynminer.py
tag=latest

ns=metadynminer-ns

build:
	docker build -t ${image}:${tag} .
	docker push ${image}:${tag}

restart: uninstall install 


install:
	helm install metadynminer -n ${ns} jupyterhub/jupyterhub \
		-f helm/values.yaml \
		--set hub.config.GenericOAuthenticator.client_secret=${shell cat helm/client_secret} \
		--set-file hub.extraConfig.form-0=helm/form-0.py \
		--set-file hub.extraConfig.pre-spawn-hook=helm/pre-spawn-hook.py \
		--set hub.config.notebookImage=${image}:${tag} \


uninstall:
	helm uninstall metadynminer -n ${ns}

log:
	kubectl -n ${ns} logs -f $(shell kubectl -n metadynminer-ns get pods | grep hub | cut -f1 -d' ')


usage:
	kubectl -n ${ns} exec $(shell kubectl -n metadynminer-ns get pods | grep hub | cut -f1 -d' ') -- cat /srv/jupyterhub/usage.log
