#!/usr/bin/env python3

has_qt6 = False
try:
    from PyQt6.QtWidgets import QApplication
    has_qt6 = True
except ModuleNotFoundError:
    print("Qt6 not found. Using Qt5 fallback")
    from PyQt5.QtWidgets import QApplication
from GUI import GUI
from Constants import Constants
from Settings import Settings
from compute.backend import LocalBackend
import sys

class QPMsimitar:
    def __init__(self):
        self.config = Settings()
        self.config.standardSettings()
        self.config.loadSettings()
        constants=Constants()
        print(constants.pi)

        # Set up compute backend based on config
        backend_type = self.config.get("Compute Backend", "local")
        if backend_type == "remote":
            try:
                from compute.backend import RemoteBackend
                self.backend = RemoteBackend(
                    server_url=self.config.get("Remote Server URL"),
                    api_token=self.config.get("Remote API Token"),
                    client_cert=self.config.get("Remote Client Cert"),
                    client_key=self.config.get("Remote Client Key"),
                    ca_cert=self.config.get("Remote CA Cert"),
                    verify_ssl=self.config.get("Remote Verify SSL") is not False
                )
            except ImportError:
                print("Remote backend dependencies not available, falling back to local")
                self.backend = LocalBackend()
        else:
            self.backend = LocalBackend()

        self.gui = GUI(self.config, self.backend)

    def showGUI(self):
        self.gui.showWindow()

if __name__ == '__main__':
    app = QApplication(sys.argv)
    MPA = QPMsimitar()
    if has_qt6:
        app.exec()
    else:
        app.exec_()
    MPA.config.saveSettings()
