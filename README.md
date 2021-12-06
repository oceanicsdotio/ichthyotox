## Drawing conclusions from drunk fish in dynamic environments.

### Post-processing raw data

Files are organized by particle types, `fish` and `cyanobacteria`. The suffixes `_ini.dat` and `_var.dat` are initial conditions for the model, and aren't relevant here. Location data is in `_position.dat` files, and state variables are `_state.dat`.

Environmental data are stored in `dissolved_toxin.dat`.

We'll ingest these in the `ichthytox` database


### Prepare

Create a table with a `postgis` geometry column.

```sql
CREATE TYPE particle_type as enum('fish', 'cyanobacteria');
CREATE TABLE locations(
        simulation           INTEGER NOT NULL,
        particle             particle_type,
        id                   INTEGER NOT NULL,
        time                 REAL NOT NULL,
    geometry             GEOMETRY NOT NULL,
        UNIQUE (simulation, particle, id, time)
);
CREATE TABLE simulation(
    uuid      UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
    id        INTEGER,
    dti       REAL NOT NULL DEFAULT 0.02,
    instp     REAL NOT NULL DEFAULT 1.00,
    dtout     REAL NOT NULL DEFAULT 0.10,
    dhor      REAL NOT NULL DEFAULT 0.10,
    dtrw      REAL NOT NULL DEFAULT 0.02,
    tdrift    INTEGER NOT NULL,
    yearlag   INTEGER NOT NULL,
    monthlag  INTEGER NOT NULL,
    daylag    INTEGER NOT NULL,
    hourlag   INTEGER NOT NULL DEFAULT 0,
    irw       INTEGER NOT NULL DEFAULT 0,
    p_sigma   BOOL NOT NULL DEFAULT false,
    out_sigma BOOL NOT NULL DEFAULT false,
    f_depth   BOOL NOT NULL DEFAULT false,
    geoarea   VARCHAR(50) NOT NULL,
    infofile  VARCHAR(50) NOT NULL DEFAULT 'screen',
    inpdir    VARCHAR(50) NOT NULL DEFAULT '/',
    lagini    VARCHAR(50) NOT NULL DEFAULT '/',
    outdir    VARCHAR(50) NOT NULL DEFAULT '/'
);
```

```
CREATE TABLE features_of_interest(
    name     VARCHAR(50) NOT NULL,
        uuid     UUID PRIMARY KEY DEFAULT uuid_generate_v4(),
        geometry INTEGER ARRAY
);
CREATE TABLE datastreams(
        things          UUID NOT NULL references things(uuid),
        time            REAL NOT NULL,
    carbohydrate    REAL,
    protein         REAL,
        mass            REAL,
    microcystin     REAL,
    UNIQUE (things, time)
);
```



### Ingest
```bash
split -l 1000000 particles.csv particles-csv-;
parts=$(ls particles-csv-*);
props="simulation, particle, time, id, geometry"
for part in ${parts}; do
    PGPASSWORD=${PG_PASSWORD} psql \
    --echo-queries \
    --username=${PG_USERNAME} \
    --host=${PG_HOST} \
    --port=${PG_PORT} \
    --dbname=ichthyotox \
    --command="\copy locations (${props}) FROM '${part}' WITH CSV"
    rm ${part}
done
```

### Create a compound index

```sql
CREATE INDEX idx_locations ON locations USING GIST(geometry);
```

For 34 million points, this takes about 24 minutes.


