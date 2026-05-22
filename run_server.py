import os
import ssl
import uvicorn
from compute.security import generate_mtls_certs

if __name__ == "__main__":
    certs_dir = os.path.join(os.path.dirname(__file__), "certs")
    ca_crt, server_crt, server_key, client_crt, client_key = generate_mtls_certs(certs_dir)
    
    print("\n--- mTLS Certificates Generated ---")
    print(f"Server Cert: {server_crt}")
    print(f"Client Cert: {client_crt}")
    print("Use these in your QPMsimitar Settings for Remote execution.")
    print("-----------------------------------\n")

    # Set default API token for testing if not set
    if "QPMSIMITAR_API_TOKEN" not in os.environ:
        os.environ["QPMSIMITAR_API_TOKEN"] = "dev-token-1234"
        print("Using default API token: 'dev-token-1234'")

    # Configure uvicorn for mTLS
    ssl_context = ssl.create_default_context(ssl.Purpose.CLIENT_AUTH)
    ssl_context.load_cert_chain(certfile=str(server_crt), keyfile=str(server_key))
    ssl_context.load_verify_locations(cafile=str(ca_crt))
    ssl_context.verify_mode = ssl.CERT_REQUIRED

    print("Starting QPMsimitar Compute Server on https://localhost:8443...")
    
    uvicorn.run("compute.server:app", host="127.0.0.1", port=8443, 
                ssl_keyfile=str(server_key), 
                ssl_certfile=str(server_crt),
                ssl_ca_certs=str(ca_crt),
                ssl_cert_reqs=ssl.CERT_REQUIRED)
