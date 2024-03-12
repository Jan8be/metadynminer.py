image=cerit.io/ljocha/metadynminer
tag=latest

ns=krenek-ns

build:
	docker build -t ${image}:${tag} .
	docker push ${image}:${tag}



install:
	helm install metadynminer -n ${ns} jupyterhub/jupyterhub \
		-f helm/values.yaml \
		--set-file hub.extraConfig.form-0=helm/form-0.py \
		--set-file hub.extraConfig.pre-spawn-hook=helm/pre-spawn-hook.py \
		--set hub.config.notebookImage=${image}:${tag} \


uninstall:
	helm uninstall metadynminer -n ${ns}
