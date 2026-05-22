import os
import secrets
import subprocess
from pathlib import Path
from fastapi import Security, HTTPException, status
from fastapi.security import APIKeyHeader

API_KEY_NAME = "X-API-Token"
api_key_header = APIKeyHeader(name=API_KEY_NAME, auto_error=False)

def get_api_token() -> str:
    """Get the expected API token from environment variable."""
    # This should be set on the server before starting
    return os.environ.get("QPMSIMITAR_API_TOKEN", "default_secret_token")

async def verify_token(api_key_header: str = Security(api_key_header)):
    """FastAPI dependency to verify the API token."""
    expected = get_api_token()
    if not api_key_header or not secrets.compare_digest(api_key_header, expected):
        raise HTTPException(
            status_code=status.HTTP_403_FORBIDDEN, detail="Invalid or missing API token"
        )
    return api_key_header

def generate_mtls_certs(out_dir: str, san_list: list = None):
    """
    Generate CA, Server, and Client certificates for mTLS using openssl.
    Allows specifying a list of Subject Alternative Names (SANs) for the server.
    Optimized to preserve CA and Client certificates while always regenerating
    the Server certificate on startup to match current host IPs.
    """
    out = Path(out_dir)
    out.mkdir(parents=True, exist_ok=True)
    
    ca_key = out / "ca.key"
    ca_crt = out / "ca.crt"
    server_key = out / "server.key"
    server_csr = out / "server.csr"
    server_crt = out / "server.crt"
    client_key = out / "client.key"
    client_csr = out / "client.csr"
    client_crt = out / "client.crt"
    
    # Compile SAN string
    default_sans = ["DNS:localhost", "IP:127.0.0.1"]
    if san_list:
        for s in san_list:
            if s and s not in default_sans:
                default_sans.append(s)
    san_str = ", ".join(default_sans)

    # Temporary config files for signing extensions
    ca_ext_file = out / "ca_ext.cnf"
    server_ext_file = out / "server_ext.cnf"
    client_ext_file = out / "client_ext.cnf"
    
    ca_ext_file.write_text("[v3_ca]\nbasicConstraints = critical, CA:true\nkeyUsage = critical, keyCertSign, cRLSign\n")
    server_ext_file.write_text(f"[v3_server]\nbasicConstraints = CA:false\nkeyUsage = critical, digitalSignature, keyEncipherment\nextendedKeyUsage = serverAuth\nsubjectAltName = {san_str}\n")
    client_ext_file.write_text("[v3_client]\nbasicConstraints = CA:false\nkeyUsage = critical, digitalSignature\nextendedKeyUsage = clientAuth\n")
    
    # 1. Generate CA (only if missing)
    if not (ca_key.exists() and ca_crt.exists()):
        print(f"Generating new Certificate Authority (CA) in {out_dir}...")
        subprocess.run(["openssl", "req", "-x509", "-newkey", "rsa:4096", "-days", "3650", 
                        "-nodes", "-keyout", str(ca_key), "-out", str(ca_crt), 
                        "-subj", "/CN=QPMsimitarCA", "-extensions", "v3_ca", "-config", str(ca_ext_file)], check=True)
                    
    # 2. Generate Client Cert (only if missing)
    if not (client_key.exists() and client_crt.exists()):
        print(f"Generating new Client Certificate in {out_dir}...")
        subprocess.run(["openssl", "req", "-newkey", "rsa:4096", "-nodes", 
                        "-keyout", str(client_key), "-out", str(client_csr), 
                        "-subj", "/CN=QPMsimitarClient"], check=True)
        subprocess.run(["openssl", "x509", "-req", "-in", str(client_csr), "-CA", str(ca_crt), 
                        "-CAkey", str(ca_key), "-CAcreateserial", "-out", str(client_crt), 
                        "-days", "365", "-sha256", "-extfile", str(client_ext_file), "-extensions", "v3_client"], check=True)

    # 3. Generate Server Cert (ALWAYS regenerate to match the current IP/SAN list)
    print(f"Regenerating Server Certificate with SANs [{san_str}]...")
    # Delete old server files if they exist to prevent prompt blocks
    for f in (server_key, server_csr, server_crt):
        try:
            f.unlink()
        except OSError:
            pass

    subprocess.run(["openssl", "req", "-newkey", "rsa:4096", "-nodes", 
                    "-keyout", str(server_key), "-out", str(server_csr), 
                    "-subj", "/CN=localhost"], check=True)
    subprocess.run(["openssl", "x509", "-req", "-in", str(server_csr), "-CA", str(ca_crt), 
                    "-CAkey", str(ca_key), "-CAcreateserial", "-out", str(server_crt), 
                    "-days", "365", "-sha256", "-extfile", str(server_ext_file), "-extensions", "v3_server"], check=True)
                    
    # Clean up temp configuration files and CA serial
    for f in (ca_ext_file, server_ext_file, client_ext_file, out / "ca.srl"):
        try:
            f.unlink()
        except OSError:
            pass
            
    return ca_crt, server_crt, server_key, client_crt, client_key
