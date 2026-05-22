import os
import ssl
import uvicorn
import secrets
import argparse
from compute.security import generate_mtls_certs

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Start secure QPMsimitar remote compute server.")
    parser.add_argument("--host", default="127.0.0.1", help="Host IP to bind the server to (default: 127.0.0.1)")
    parser.add_argument("--port", type=int, default=8443, help="Port to listen on (default: 8443)")
    parser.add_argument("--san", action="append", help="Additional Subject Alternative Names (DNS or IP) for the server certificate. Can be specified multiple times.")
    args = parser.parse_args()

    import socket
    san_list = []
    
    # Auto-detect host type and add to SAN list
    host = args.host
    if host not in ("127.0.0.1", "localhost", "0.0.0.0"):
        try:
            socket.inet_aton(host)
            san_list.append(f"IP:{host}")
        except socket.error:
            san_list.append(f"DNS:{host}")
            
    # Add any user specified SANs
    if args.san:
        for s in args.san:
            if not (s.startswith("IP:") or s.startswith("DNS:")):
                try:
                    socket.inet_aton(s)
                    san_list.append(f"IP:{s}")
                except socket.error:
                    san_list.append(f"DNS:{s}")
            else:
                san_list.append(s)

    certs_dir = os.path.join(os.path.dirname(__file__), "certs")
    ca_crt, server_crt, server_key, client_crt, client_key = generate_mtls_certs(certs_dir, san_list=san_list)
    
    print("\n--- mTLS Certificates Generated ---")
    print(f"Server Cert: {server_crt}")
    print(f"Client Cert: {client_crt}")
    print("Use these in your QPMsimitar Settings for Remote execution.")
    print("-----------------------------------\n")

    # Set default random API token if not set in environment
    if "QPMSIMITAR_API_TOKEN" not in os.environ:
        random_token = secrets.token_hex(16)
        os.environ["QPMSIMITAR_API_TOKEN"] = random_token
        print(f"Generated secure, random API token for this session:\n>>> {random_token} <<<\n")
    else:
        print(f"Using pre-configured API token from environment: '{os.environ['QPMSIMITAR_API_TOKEN']}'")

    # Configure uvicorn for mTLS
    ssl_context = ssl.create_default_context(ssl.Purpose.CLIENT_AUTH)
    ssl_context.load_cert_chain(certfile=str(server_crt), keyfile=str(server_key))
    ssl_context.load_verify_locations(cafile=str(ca_crt))
    ssl_context.verify_mode = ssl.CERT_REQUIRED

    print(f"Starting QPMsimitar Compute Server on https://{args.host}:{args.port}...")
    
    uvicorn.run("compute.server:app", host=args.host, port=args.port, 
                ssl_keyfile=str(server_key), 
                ssl_certfile=str(server_crt),
                ssl_ca_certs=str(ca_crt),
                ssl_cert_reqs=ssl.CERT_REQUIRED)
