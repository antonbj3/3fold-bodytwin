#!/usr/bin/env bash
# Stop BodyTwin's queue driver (the bt-queue-driver service in jobs-bodytwin.slice); leave started jobs untouched.
systemctl --user stop bt-queue-driver.service 2>/dev/null && echo "bt-queue-driver stopped"
