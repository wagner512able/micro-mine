# Guia rápido (português): GitHub + Vercel

1. **Teste no seu computador**: `pip install -e ".[web]"` e depois `python scripts/dev_vercel.py`; abra http://127.0.0.1:3000.
2. **GitHub**: crie uma conta e um repositório público chamado `micro-mine`. Envie os arquivos
   (pelos comandos `git` do arquivo `DEPLOY_VERCEL.md` ou arrastando a pasta em *Add file > Upload files*).
   Troque `YOUR-USERNAME` pelo seu usuário nos links.
3. **Vercel**: entre em vercel.com com a conta do GitHub, clique em *Add New > Project*, escolha o repositório,
   deixe *Framework Preset = Other* e clique em *Deploy*. Em alguns minutos o site recebe um endereço
   `https://micro-mine-xxxx.vercel.app`.
4. **Atualizações**: cada `git push` na branch `main` publica de novo sozinho. Se alterar a página, rode antes
   `python scripts/sync_public.py`.
5. **Limites da versão gratuita online**: cerca de 4 MB de FASTA comprimido por análise (em torno de 50 cloroplastos de
   155 kb) e 60 s de processamento. Para jobs maiores, use o `Dockerfile` em um servidor, ou instale localmente.
6. **Licença**: o arquivo `LICENSE` usa MIT como sugestão; confirme a política da sua instituição antes de publicar.
