# Steering LANE_FOCAL_DELIVERY_GEOMETRY — new lane 2026-10-04

## Starting point
The decision already exists as a runnable file: `tasks/assembly/focal_delivery_decision.py`.
Run it first. Read what it does. Build on it, do not redo it.

## What it says
Two independent measurement pairs give the same diffusivity within 1,40 times, 0,175 to 0,245 µm² per second. A 10-second pulse carries 3,13 µm. That is the minimum allowed distance between delivery points.

## The obstacle
The decision has been made but not tested against the world. No measurement has confirmed or rejected it.

## Your task
Find the measured ground truth that decides it. Published, with DOI or PMID.
If ground truth exists: recalculate our number against it and say whether the decision stands.
If no ground truth exists: specify the measurement. Quantity, unit, resolution, number of samples.
A specification of what is missing is a result.

## Control
The equally informed reading that is in the file. Read it before building your own.
The gain counts only against it.

## Falsifier
Below that distance, point delivery and a uniform bath are the same stimulus. Activation ratio ≤ 1 above the boundary rejects the decision.

## Rules
Unit in every field name. Plane, convention and measure of spread written out.
Six collisions in the project in one day came from names that did not carry their convention.
Everything PENDING_INDEPENDENT_REVIEW. No breakthroughs without an equally informed control.
